// @(#)root/utils:$Id$
// Author: Philippe Canal 27/08/2003

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun, Fons Rademakers, and al.          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS              *
 *************************************************************************/

#include "RConfigure.h"
#include <ROOT/RConfig.hxx>

#include "RStl.h"
#include "TClassEdit.h"
#include "TClingUtils.h"
using namespace TClassEdit;

#include <stdio.h>

#include "clang/AST/Decl.h"
#include "clang/AST/DeclTemplate.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"

#include "clang/Sema/Sema.h"
#include "clang/Sema/Template.h"
#include "clang/Frontend/CompilerInstance.h"

#include "Varargs.h"

//
// ROOT::Internal::RStl is the rootcint STL handling class.
//

static int fgCount = 0;

ROOT::Internal::RStl& ROOT::Internal::RStl::Instance()
{
   // Return the singleton ROOT::Internal::RStl.

   static ROOT::Internal::RStl instance;
   return instance;

}

void ROOT::Internal::RStl::GenerateTClassFor(const clang::QualType &type, const cling::Interpreter &interp, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
{
   // Force the generation of the TClass for the given class.

   clang::QualType thisType = type;

   auto typePtr = thisType.getTypePtr();
   const clang::CXXRecordDecl *stlclass = typePtr->getAsCXXRecordDecl();
   if (!stlclass) {
      return;
   }

   // Transform the type to the corresponding one for IO
   auto typeForIO = ROOT::TMetaUtils::GetTypeForIO(thisType, interp, normCtxt);
   if (typeForIO.getTypePtr() != typePtr)
      stlclass = typeForIO->getAsCXXRecordDecl();
   if (!stlclass) {
      return;
   }
   thisType = typeForIO;

   const clang::ClassTemplateSpecializationDecl *templateCl = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(stlclass);

   if (!templateCl) {
      ROOT::TMetaUtils::Error("RStl::GenerateTClassFor","%s not in a template",
            ROOT::TMetaUtils::GetQualifiedName(*stlclass).c_str());
   }

   if ( TClassEdit::STLKind( stlclass->getName().str() )  == ROOT::kSTLvector ) {
      const clang::TemplateArgument &arg( templateCl->getTemplateArgs().get(0) );
      if (arg.getKind() == clang::TemplateArgument::Type) {
         const clang::NamedDecl *decl = arg.getAsType().getTypePtr()->getAsCXXRecordDecl();
         if (decl) {
            // NOTE: we should just compare the decl to the bool builtin!
            llvm::StringRef argname = decl->getName();
            if ( (argname.str() == "bool") || (argname.str() == "Bool_t") ) {
               ROOT::TMetaUtils::Warning("std::vector<bool>", " is not fully supported yet!\nUse std::vector<char> or std::deque<bool> instead.\n");
            }
         }
      }
   }

   fList.insert( ROOT::TMetaUtils::AnnotatedRecordDecl(++fgCount,
                                                       thisType.getTypePtr(),
                                                       stlclass,
                                                       "",
                                                       false /* for backward compatibility rather than 'true' .. neither really make a difference */,
                                                       false,
                                                       false,
                                                       false,
                                                       -1,
                                                       interp,
                                                       normCtxt) );

   // fprintf(stderr,"Registered the STL class %s as needing a dictionary\n",R__GetQualifiedName(*stlclass).c_str());

   for(unsigned int i=0; i <  templateCl->getTemplateArgs().size(); ++i) {
      const clang::TemplateArgument &arg( templateCl->getTemplateArgs().get(i) );
      if (arg.getKind() == clang::TemplateArgument::Type) {
         const clang::NamedDecl *decl = arg.getAsType().getTypePtr()->getAsCXXRecordDecl();

         if (decl && TClassEdit::STLKind( decl->getName().str() ) != 0 )
            {
               const clang::CXXRecordDecl *clxx = llvm::dyn_cast<clang::CXXRecordDecl>(decl);
               if (clxx) {
                  if (!clxx->isCompleteDefinition()) {
                     /* bool result = */ ROOT::TMetaUtils::RequireCompleteType(interp, clxx->getLocation(), arg.getAsType());
                  }
                  // Do we need to strip the qualifier?
                  GenerateTClassFor(arg.getAsType(),interp,normCtxt);
               }
            }
      }
   }
}

void ROOT::Internal::RStl::GenerateTClassFor(const char *requestedName, const clang::CXXRecordDecl *stlclass, const cling::Interpreter &interp, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
{
   // Force the generation of the TClass for the given class.
   const clang::ClassTemplateSpecializationDecl *templateCl = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(stlclass);

   if (templateCl == 0) {
      ROOT::TMetaUtils::Error("RStl::GenerateTClassFor","%s not in a template",
            ROOT::TMetaUtils::GetQualifiedName(*stlclass).c_str());
   }


   if ( TClassEdit::STLKind( stlclass->getName().str() )  == ROOT::kSTLvector ) {
      const clang::TemplateArgument &arg( templateCl->getTemplateArgs().get(0) );
      if (arg.getKind() == clang::TemplateArgument::Type) {
         const clang::NamedDecl *decl = arg.getAsType().getTypePtr()->getAsCXXRecordDecl();
         if (decl) {
            // NOTE: we should just compare the decl to the bool builtin!
            llvm::StringRef argname = decl->getName();
            if ( (argname.str() == "bool") || (argname.str() == "Bool_t") ) {
               ROOT::TMetaUtils::Warning("std::vector<bool>", " is not fully supported yet!\nUse std::vector<char> or std::deque<bool> instead.\n");
            }
         }
      }
   }

   fList.insert( ROOT::TMetaUtils::AnnotatedRecordDecl(++fgCount,stlclass,requestedName,true,false,false,false,-1, interp,normCtxt) );

   TClassEdit::TSplitType splitType( requestedName, (TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd) );
   for(unsigned int i=0; i <  templateCl->getTemplateArgs().size(); ++i) {
      const clang::TemplateArgument &arg( templateCl->getTemplateArgs().get(i) );
      if (arg.getKind() == clang::TemplateArgument::Type) {
         const clang::NamedDecl *decl = arg.getAsType().getTypePtr()->getAsCXXRecordDecl();

         if (decl && TClassEdit::STLKind( decl->getName().str() ) != 0 )
            {
               const clang::CXXRecordDecl *clxx = llvm::dyn_cast<clang::CXXRecordDecl>(decl);
               if (clxx) {
                  if (!clxx->isCompleteDefinition()) {
                     /* bool result = */ ROOT::TMetaUtils::RequireCompleteType(interp, clxx->getLocation (), arg.getAsType());
                     clxx = arg.getAsType().getTypePtr()->getAsCXXRecordDecl();
                  }
                  if (!splitType.fElements.empty()) {
                     GenerateTClassFor( splitType.fElements[i+1].c_str(), clxx, interp, normCtxt);
                  } else {
                     GenerateTClassFor( "", clxx, interp, normCtxt );
                  }
               }
            }
      }
   }

}

void ROOT::Internal::RStl::Print()
{
   // Print the content of the object
   fprintf(stderr,"ROOT::Internal::RStl singleton\n");
   list_t::iterator iter;
   for(iter = fList.begin(); iter != fList.end(); ++iter) {
      fprintf(stderr, "need TClass for %s\n", ROOT::TMetaUtils::GetQualifiedName(*(*iter)).c_str());
   }
}

void ROOT::Internal::RStl::WriteClassInit(std::ostream &ostr,
                                const cling::Interpreter &interp,
                                const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt,
                                const ROOT::TMetaUtils::RConstructorTypes& ctorTypes,
                                bool &needCollectionProxy,
                                void (*emitStreamerInfo)(const char*) )
{
   // This function writes the TGeneraticClassInfo initialiser
   // and the auxiliary functions (new and delete wrappers) for
   // each of the STL containers that have been registered

   list_t::iterator iter;
   for(iter = fList.begin(); iter != fList.end(); ++iter) {
      const clang::CXXRecordDecl* result;

      if (!iter->GetRecordDecl()->getDefinition()) {

         // We do not have a complete definition, we need to force the instantiation
         // and findScope can do that.
         const cling::LookupHelper& lh = interp.getLookupHelper();
         result = llvm::dyn_cast_or_null<clang::CXXRecordDecl>(lh.findScope(iter->GetNormalizedName(),
                                                                            cling::LookupHelper::NoDiagnostics,
                                                                            0)
                                                               );

         if (!result || !iter->GetRecordDecl()->getDefinition()) {
            fprintf(stderr,"Error: incomplete definition for %s\n",iter->GetNormalizedName());
            continue;
         }
      }
      else
      {
         result = llvm::dyn_cast<clang::CXXRecordDecl>(iter->GetRecordDecl());
      }

      ROOT::TMetaUtils::WriteClassInit(ostr, *iter, result, interp, normCtxt, ctorTypes, needCollectionProxy);
      ROOT::TMetaUtils::WriteAuxFunctions(ostr, *iter, result, interp, ctorTypes, normCtxt);

      if (emitStreamerInfo) emitStreamerInfo(iter->GetNormalizedName());
   }
}
