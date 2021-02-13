// @(#)root/utils/src:$Id$
// Author: Danilo Piparo January 2014

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/rootcint.            *
 *************************************************************************/

#ifndef __DICTSELECTIONREADER__
#define __DICTSELECTIONREADER__

#include "clang/AST/RecursiveASTVisitor.h"

#include <llvm/ADT/StringMap.h>

#include <set>
#include <unordered_set>
#include <string>
#include <unordered_map>

class SelectionRules;
class ClassSelectionRule;
namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}
namespace cling {
class Interpreter;
}

namespace clang {
   class ASTContext;
//    class DeclContext;
   class NamespaceDecl;
   class CXXRecordDecl;
}

/**
 * @file  DictSelectionReader.h
 * @author Danilo Piparo
 * @date January 2014
 * @brief Select classes and assign properties using C++ syntax.
 *
 * @details When generating dictionary information for a class,
 * one sometimes wants to specify additional information
 * beyond the class definition itself, for example, to specify
 * that certain members are to be treated as transient by the persistency
 * system.  This can be done by associating a dictionary selection class
 * with the class for which dictionary information is being generated.
 * The contents of this selection class encode the additional information.
 * Below, we first discuss how to associate a selection class
 * with your class; then we list the current Set of information
 * which may appear inside the selection class.
 *
 * The simplest case is for the case of a non-template class @c C.
 * By default, the Name of the selection class is then
 * @c ROOT::Meta::Selection::C.  If you have such a class, it will be found
 * automatically.  If @c C is in a namespace, @c NS::C, then
 * the selection class should be in the same namespace:
 * @c ROOT::Meta::Selection::NS::C.
 *
 * Examples:
 *
 * The DictSelectionReader is used to create selection rules starting from
 * C++ the constructs of the @c ROOT::Meta::Selection namespace. All rules
 * are matching by name.
 * A brief description of the operations that lead to class selection:
 * 1. If a class declaration is present in the selection namespace, a class
 *    with the same name is selected outside the selection namespace.
 * 2. If a template class declaration and a template instantiation is present
 *    in the selection namespace, all the instances of the template are
 *    selected outside the namespace.
 *
 * For example:
 * @code
 * [...]
 * class classVanilla{};
 * template <class A> class classTemplateVanilla {};
 * classTemplateVanilla<char> t0;
 * namespace ROOT {
 *    namespace Meta {
 *       namespace Selection {
 *          class classVanilla {};
 *          template <typename A> class classTemplateVanilla {};
 *          classTemplateVanilla<char> st0;
 *       }
 *    }
 * }
 * @endcode
 * would create two selection rules to select @c classVanilla and
 * @c classTemplateVanilla<char>.
 *
 * A brief description of the properties that can be assigned to classes
 * with the @c ROOT::Meta::Selection::ClassAttributes class.
 * 1. @c kNonSplittable : Makes the class non splittable
 *
 * The class properties can be assigned via a traits mechanism. For example:
 * @code
 * [...]
 * class classWithAttributes {};
 * namespace ROOT {
 *    namespace Meta {
 *       namespace Selection {
 *          class classWithAttributes : ClassAttributes <kNonSplittable> {};
 *       }
 *    }
 * }
 * @endcode
 * would create a selection rule which selects class @c classWithAttributes and
 * assignes to it the property described by @c kNonSplittable. Multiple
 * properties can be assigned to a single class with this syntax:
 * @code
 * [...]
 * namespace ROOT {
 *    namespace Meta {
 *       namespace Selection {
 *          class classWithAttributes :
 *             ClassAttributes <kProperty1 + kProperty2 + ... + kPropertyN> {};
 *       }
 *    }
 * }
 * @endcode
 *
 * The @c ROOT::Meta::Selection syntax allows to alter the number of template
 * parameters of a certain template class within the ROOT type system, TClass.
 * Technically it allows to alter the way in which the "normalized name" (in
 * other words, the "ROOT name") of the class is created. The key is the usage
 * of the @c KeepFirstTemplateArguments traits class.
 * It is possible to select the maximum number of template arguments considered
 * if not different from the default. A concrete example can be more clear than
 * a long explaination in this case:
 * @code
 * [...]
 * template <class T, class U = int, int V = 3> class A {...};
 * template <class T, class Alloc = myAllocator<T> > class myVector {...};
 * A<char> a1;
 * A<char, float> a2;
 * myVector<float> v1;
 * myVector<A<char>> v2;
 *
 * namespace ROOT {
 *    namespace Meta {
 *       namespace Selection {
 *          template <class T, class U = int, int V = 3> class A
 *            : KeepFirstTemplateArguments<1> {};
 *          A<double> a;
 *          template <class T, class Alloc = myAllocator<T> > class myVector
 *            : KeepFirstTemplateArguments<1> {};
 *          myVector<double> vd;
 *       }
 *    }
 * }
 * @endcode
 *
 * Consistently with what described above, all the instances of @c A and
 * @c myvector will be selected. In addition, only the first template parameter
 * will be kept.
 * In absence of any @c KeepFirstTemplateArguments trait, the normalization
 * would be:
 * - @c A<char>           &rarr; @c A<char,float,3>
 * - @c A<char,float>     &rarr; @c A<char,int,3>
 * - @c myVector<float>   &rarr; @c myVector<A<char,int,3>,myAllocator<A<char,int,3>>>
 * - @c myVector<A<char>> &rarr; @c myVector<float,myAllocator<float>>
 *
 * Now, deciding to keep just one argument (@c KeepFirstTemplateArguments<1>):
 * - @c A<char>           &rarr; @c A<char,float>
 * - @c A<char,float>     &rarr; @c A<char>
 * - @c myVector<float>   &rarr; @c myVector<A<char>,myAllocator<A<char>>>
 * - @c myVector<A<char>> &rarr; @c myVector<float,myAllocator<float>>
 *
 * And deciding to keep two arguments (@c KeepFirstTemplateArguments<2>):
 * - @c A<char>           &rarr; @c A<char,float>
 * - @c A<char,float>     &rarr; @c A<char,int>
 * - @c myVector<float>   &rarr; @c myVector<A<char,int>,myAllocator<A<char,int>>>
 * - @c myVector<A<char>> &rarr; @c myVector<float,myAllocator<float>>
 *
 * A brief description of the properties that can be assigned to data members
 * with the @c ROOT::Meta::Selection MemberAttributes class:
 * 1. @c kTransient : the data member is transient, not persistified by the
 *    ROOT I/O.
 * 2. @c kAutoSelected : the type of the data member is selected without the
 *    need of specifying its class explicitely.
 *
 * For example:
 * @code
 * [...]
 * class classTransientMember {
 *  private:
 *    int transientMember;
 * };
 * class classAutoselected {};
 * class classTestAutoselect {
 *  private:
 *    classAutoselected autoselected;
 * };
 *
 * namespace ROOT {
 *    namespace Meta {
 *       namespace Selection {
 *          class classTestAutoselect {
 *             MemberAttributes<kAutoSelected> autoselected;
 *          };
 *          class classTransientMember {
 *             MemberAttributes<kTransient> transientMember;
 *          };
 *       }
 *    }
 * }
 * @endcode
 * would lead to the creation of selection rules for @c classTransientMember
 * specifying that @c transientMember is transient, @c classTestAutoselect and
 * @c classAutoselected.
 *
 * Another trait class present in the @c ROOT::Meta::Selection is
 * @c SelectNoInstance. If a template in the selection namespace inherits from
 * this class, none of its instantiations will be automatically selected but
 * all of the properties specified otherwise, like transient members or
 * number of template arguments to keep, will be transmitted to all of the
 * instantiations selected by other means.
 * For example
 * @code
 * [...]
 * template< class T, class BASE >
 * class MyDataVector : KeepFirstTemplateArguments< 1 >, SelectNoInstance {
 *     MemberAttributes< kTransient + kAutoSelected > m_isMostDerived;
 *     MemberAttributes< kNonSplittable + kAutoSelected > m_isNonSplit;
 * };
 * [...]
 * @endcode
 **/
namespace ROOT {
namespace Internal {

class DictSelectionReader : public clang::RecursiveASTVisitor<DictSelectionReader> {
public:
   /// Take the selection rules as input (for consistency w/ other selector
   /// interfaces)
   DictSelectionReader(cling::Interpreter &interp, SelectionRules &, const clang::ASTContext &,
                       ROOT::TMetaUtils::TNormalizedCtxt &);

   /// Visit the entities that needs to be selected
   bool VisitRecordDecl(clang::RecordDecl *);

   bool shouldVisitTemplateInstantiations() const {
      return true;
   }

private:


   struct TemplateInfo { /// < Class to store the information about templates upon parsing
      TemplateInfo(int argsToKeep): fArgsToKeep(argsToKeep) {};
      TemplateInfo() {};
      int fArgsToKeep = -1;
      std::unordered_set<std::string> fTransientMembers {};
      std::unordered_set<std::string> fUnsplittableMembers {};
   };

   inline bool
   InSelectionNamespace(const clang::RecordDecl &,
                        const std::string &str =
                           ""); ///< Check if in the ROOT::Meta::Selection namespace
   inline bool FirstPass(const clang::RecordDecl &); ///< First pass on the AST
   inline bool SecondPass(const clang::RecordDecl &); ///< Second pass on the AST, using the information of the first one
   inline void
   ManageFields(const clang::RecordDecl &,
                const std::string &,
                ClassSelectionRule &,
                bool); ///< Take care of the class fields
   inline void
   ManageBaseClasses(const clang::CXXRecordDecl &, const std::string &, bool &); ///< Take care of the class bases
   template <class T>
   inline unsigned int ExtractTemplateArgValue(
      const T &,
      const std::string &); ///< Extract the value of the template parameter
   inline const clang::TemplateArgumentList *GetTmplArgList(
      const clang::CXXRecordDecl &); ///< Get the template arguments list if any

   std::string PatternifyName(const std::string &className); ///< Transform instance name in pattern for selection
   void GetPointeeType(std::string &typeName); ///< Get name of the pointee type

   SelectionRules &fSelectionRules; ///< The selection rules to be filled
   std::set<const clang::RecordDecl *>
   fSelectedRecordDecls; ///< The pointers of the selected RecordDecls
   std::set<std::string>
   fSpecialNames; ///< The names of the classes used for the selction syntax
   llvm::StringMap<std::set<std::string> >
   fAutoSelectedClassFieldNames; ///< Collect the autoselected classes
   llvm::StringMap<std::set<std::string> >
   fNoAutoSelectedClassFieldNames; ///< Collect the autoexcluded classes
   std::unordered_map<std::string, TemplateInfo> fTemplateInfoMap; ///< List template name - properties map
   llvm::StringMap<ClassSelectionRule>
   fClassNameSelectionRuleMap; ///< Map of the already built sel rules
   bool fIsFirstPass; ///< Keep trance of the number of passes through the AST
   ROOT::TMetaUtils::TNormalizedCtxt &fNormCtxt; ///< The reference to the normalized context
};
}
}

#endif
