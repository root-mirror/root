#include "DictSelectionReader.h"

#include "clang/AST/AST.h"

#include "cling/Interpreter/Interpreter.h"

#include "ClassSelectionRule.h"
#include "SelectionRules.h"
#include "TClingUtils.h"
#include "TClassEdit.h"

#include "RootMetaSelection.h"

#include <iostream>
#include <sstream>

namespace ROOT {
namespace Internal {

////////////////////////////////////////////////////////////////////////////////

DictSelectionReader::DictSelectionReader(cling::Interpreter &interp, SelectionRules &selectionRules,
                                         const clang::ASTContext &C, ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
   : fSelectionRules(selectionRules), fIsFirstPass(true), fNormCtxt(normCtxt)
{
   clang::TranslationUnitDecl *translUnitDecl = C.getTranslationUnitDecl();

   {
      // We push a new transaction because we could deserialize decls here
      cling::Interpreter::PushTransactionRAII RAII(&interp);
      // Inspect the AST
      TraverseDecl(translUnitDecl);
   }

   // Now re-inspect the AST to find autoselected classes (double-tap)
   fIsFirstPass = false;
   if (!fTemplateInfoMap.empty() ||
         !fAutoSelectedClassFieldNames.empty() ||
         !fNoAutoSelectedClassFieldNames.empty())
      TraverseDecl(translUnitDecl);

   // Now push all the selection rules
   for (llvm::StringMap<ClassSelectionRule>::iterator it =
            fClassNameSelectionRuleMap.begin();
         it != fClassNameSelectionRuleMap.end();
         ++it) {
      fSelectionRules.AddClassSelectionRule(it->second);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If it's not contained by 2 namespaces, drop it.

/**
 * Check that the recordDecl is enclosed in the ROOT::Meta::Selection namespace,
 * excluding the portion dedicated the definition of the syntax, which is part
 * of ROOT, not of the user code.
 * If performance is needed, an alternative approach to string comparisons
 * could be adopted. One could use for example hashes of strings in first
 * approximation.
 **/
bool
DictSelectionReader::InSelectionNamespace(const clang::RecordDecl &recordDecl,
      const std::string &className)
{
   std::list<std::pair<std::string, bool> > enclosingNamespaces;
   ROOT::TMetaUtils::ExtractEnclosingNameSpaces(recordDecl,
         enclosingNamespaces);

   const unsigned int nNs = enclosingNamespaces.size();
   if (nNs < 3) return false;

   if (enclosingNamespaces.back().second || // is inline namespace
         enclosingNamespaces.back().first != "ROOT")
      return false;

   enclosingNamespaces.pop_back();
   if (enclosingNamespaces.back().second || // is inline namespace
         enclosingNamespaces.back().first != "Meta")
      return false;

   enclosingNamespaces.pop_back();
   if (enclosingNamespaces.back().second || // is inline namespace
         enclosingNamespaces.back().first != "Selection")
      return false;

   // Exclude the special names identifying the entities of the selection syntax
   if (className != "" &&
         (className.find("MemberAttributes") == 0 ||
          className.find("ClassAttributes") == 0 || className.find("Keep") == 0))
      return false;

   return true;
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Get the pointer to the template arguments list. Return zero if not available.
 **/
const clang::TemplateArgumentList *
DictSelectionReader::GetTmplArgList(const clang::CXXRecordDecl &cxxRcrdDecl)
{
   const clang::ClassTemplateSpecializationDecl *tmplSpecDecl =
      llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(&cxxRcrdDecl);

   if (!tmplSpecDecl) return nullptr;

   return &tmplSpecDecl->getTemplateArgs();
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Extract the value of the integral template parameter of a CXXRecordDecl when
 * it has a certain name. If nothing can be extracted, the value of @c zero
 * is returned.
 **/
template <class T>
unsigned int
DictSelectionReader::ExtractTemplateArgValue(const T &myClass,
      const std::string &pattern)
{
   const clang::RecordDecl *rcrdDecl =
      ROOT::TMetaUtils::GetUnderlyingRecordDecl(myClass.getType());
   const clang::CXXRecordDecl *cxxRcrdDecl =
      llvm::dyn_cast<clang::CXXRecordDecl>(rcrdDecl);

   if (!cxxRcrdDecl) return 0;

   const clang::TemplateArgumentList *tmplArgs = GetTmplArgList(*cxxRcrdDecl);
   if (!tmplArgs) return 0;

   if (std::string::npos == cxxRcrdDecl->getNameAsString().find(pattern))
      return 0;

   return tmplArgs->get(0).getAsIntegral().getLimitedValue();
}

////////////////////////////////////////////////////////////////////////////////
/// Iterate on the members to see if
/// 1) They are transient
/// 2) They imply further selection

/**
 * Loop over the class filelds and take actions according to their properties
 *    1. Insert a field selection rule marking a member transient
 *    2. Store in a map the name of the field the type of which should be
 * autoselected. The key is the name of the class and the value the name of the
 * field. This information is used in the second pass.
 **/
void DictSelectionReader::ManageFields(const clang::RecordDecl &recordDecl,
                                       const std::string &className,
                                       ClassSelectionRule &csr,
                                       bool autoselect)
{
   std::string pattern = className.substr(0, className.find_first_of("<"));

   for (auto fieldPtr : recordDecl.fields()) {

      unsigned int attrCode =
         ExtractTemplateArgValue(*fieldPtr, "MemberAttributes");

      if (attrCode == ROOT::Meta::Selection::kMemberNullProperty) continue;

      const char *fieldName = fieldPtr->getName().data();

      if (attrCode & ROOT::Meta::Selection::kNonSplittable) {
         if (!autoselect) {
            fTemplateInfoMap[pattern].fUnsplittableMembers.insert(fieldName);
         } else {
            VariableSelectionRule vsr(BaseSelectionRule::kYes);
            vsr.SetAttributeValue(ROOT::TMetaUtils::propNames::name, fieldName);
            vsr.SetAttributeValue(ROOT::TMetaUtils::propNames::comment, "||");
            csr.AddFieldSelectionRule(vsr);
         }
      }

      if (attrCode & ROOT::Meta::Selection::kTransient) {
         if (!autoselect) {
            fTemplateInfoMap[pattern].fTransientMembers.insert(fieldName);
         } else {
            VariableSelectionRule vsr(BaseSelectionRule::kYes);
            vsr.SetAttributeValue(ROOT::TMetaUtils::propNames::name, fieldName);
            vsr.SetAttributeValue(ROOT::TMetaUtils::propNames::comment, "!");
            csr.AddFieldSelectionRule(vsr);
         }
      }

      if (attrCode & ROOT::Meta::Selection::kAutoSelected)
         fAutoSelectedClassFieldNames[className].insert(fieldName);
      else if (attrCode & ROOT::Meta::Selection::kNoAutoSelected)
         fNoAutoSelectedClassFieldNames[className].insert(fieldName);

   } // end loop on fields
}

////////////////////////////////////////////////////////////////////////////////
/// Check the traits of the class. Useful information may be there
/// extract mothers, make a switchcase:
/// 1) templates args are to be skipped
/// 2) There are properties. Make a loop. make a switch:
///  2a) Is splittable

/**
 * Manage the loop over the base classes.
 * Initially, the class attributes are identified and selection rules filled
 * if:
 *    1. The class is not splittable
 * Then we look for the traits pointing to the need of hiding template
 * arguments. This information is stored in the form of a list of pairs, where
 * the first argument is the pattern of the template instance to match and
 * the second one the number of arguments to be skipped. This information is
 * used during the second pass.
 **/
void
DictSelectionReader::ManageBaseClasses(const clang::CXXRecordDecl &cxxRcrdDecl,
                                       const std::string &className,
                                       bool &autoselect)
{
   std::string baseName;
   clang::ASTContext &C = cxxRcrdDecl.getASTContext();
   for (auto & base : cxxRcrdDecl.bases()) {

      if (unsigned int nArgsToKeep = ExtractTemplateArgValue(base, "Keep")) {
         std::string pattern =
            className.substr(0, className.find_first_of("<"));
         // Fill the structure holding the template and the number of args to
         // skip
         fTemplateInfoMap[pattern] = TemplateInfo(nArgsToKeep);
      }

      // at most one string comparison...
      if (autoselect) {
         auto qt = base.getType();
         ROOT::TMetaUtils::GetFullyQualifiedTypeName(baseName, qt, C);
         if (baseName == "ROOT::Meta::Selection::SelectNoInstance") autoselect = false;
      }

   } // end loop on base classes
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Manage the first pass over the AST, inspecting only nodes which are within
 * the selection namespace. Selection rules are directly filled as well as
 * data sructures re-used during the second pass.
 **/
bool DictSelectionReader::FirstPass(const clang::RecordDecl &recordDecl)
{
   std::string className;
   ROOT::TMetaUtils::GetQualifiedName(
      className, *recordDecl.getTypeForDecl(), recordDecl);

   // Strip ROOT::Meta::Selection
   className.replace(0, 23, "");

   if (!InSelectionNamespace(recordDecl, className)) return true;

   if (!fSelectedRecordDecls.insert(&recordDecl).second) return true;

   bool autoselect = true;
   if (auto cxxRcrdDecl = llvm::dyn_cast<clang::CXXRecordDecl>(&recordDecl)) {
      ManageBaseClasses(*cxxRcrdDecl, className, autoselect);
   }

   ClassSelectionRule csr(BaseSelectionRule::kYes);
   const size_t lWedgePos(className.find_first_of("<"));
   std::string patternName("");
   if (lWedgePos != std::string::npos &&
         llvm::isa<clang::ClassTemplateSpecializationDecl>(recordDecl)) {
      patternName = PatternifyName(className);
      csr.SetAttributeValue(ROOT::TMetaUtils::propNames::pattern, patternName);

   } else {
      csr.SetAttributeValue(ROOT::TMetaUtils::propNames::name, className);
   }

   ManageFields(recordDecl, className, csr, autoselect);

   if (!autoselect) return true;

   // Finally add the selection rule
   fClassNameSelectionRuleMap[patternName.empty() ? className : patternName] =
      csr;

   return true;
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Second pass through the AST. Two operations are performed:
 *    1. Selection rules for classes to be autoselected are created. The
 * algorithm works as follows: the members of the classes matching the name of
 * the classes which contained autoselected members in the selection namespace
 * are inspected. If a field with the same name of the one which was
 * autoselected a selection rule based on its typename is built.
 *    2. If a class is found which is a @c TemplateSpecialisationDecl its
 * name is checked to match one of the patterns identified during the first
 * pass. If a match is found, a property is added to the selection rule with
 * the number of template arguments to keep in order to percolate this
 * information down to the @c AnnotatedRecordDecl creation which happens in the
 * @c RScanner .
 **/
bool DictSelectionReader::SecondPass(const clang::RecordDecl &recordDecl)
{
   using namespace ROOT::TMetaUtils;

   // No interest if we are in the selection namespace
   if (InSelectionNamespace(recordDecl)) return true;

   std::string className;
   GetQualifiedName(className, *recordDecl.getTypeForDecl(), recordDecl);

   // If the class is not among those which have fields the type of which are to
   // be autoselected or excluded
   if (0 != fAutoSelectedClassFieldNames.count(className) ||
         0 != fNoAutoSelectedClassFieldNames.count(className)) {
      // Iterate on fields. If the name of the field is among the ones the types
      // of which should be (no)autoselected, add a class selection rule
      std::string typeName;
      clang::ASTContext &C = recordDecl.getASTContext();
      for (clang::RecordDecl::field_iterator filedsIt =
               recordDecl.field_begin();
            filedsIt != recordDecl.field_end();
            ++filedsIt) {
         const std::string fieldName(filedsIt->getNameAsString());
         bool excluded = 1 == fNoAutoSelectedClassFieldNames[className].count(fieldName);
         bool selected = 1 == fAutoSelectedClassFieldNames[className].count(fieldName);
         if (!selected && !excluded)
            continue;
         ClassSelectionRule aSelCsr(excluded ? BaseSelectionRule::kNo : BaseSelectionRule::kYes);
         GetFullyQualifiedTypeName(typeName, filedsIt->getType(), C);
         GetPointeeType(typeName);
         aSelCsr.SetAttributeValue(propNames::name, typeName);
         fSelectionRules.AddClassSelectionRule(aSelCsr);
      }
   }

   // If the class is a template instantiation and its name matches one of the
   // patterns

   // We don't want anything different from templ specialisations
   if (auto tmplSpecDecl = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(&recordDecl)) {
      for (auto & patternInfoPair : fTemplateInfoMap) {
         const std::string &pattern = patternInfoPair.first;
         const TemplateInfo &tInfo = patternInfoPair.second;
         // Check if we have to add a selection rule for this class
         if (className.find(pattern) != 0) continue;

         // Take care of the args to keep
         auto ctd = tmplSpecDecl->getSpecializedTemplate();
         if (tInfo.fArgsToKeep != -1 && ctd) {
            fNormCtxt.AddTemplAndNargsToKeep(ctd->getCanonicalDecl(), tInfo.fArgsToKeep);
         }

         // Now we take care of the transient and unsplittable members
         if (tInfo.fTransientMembers.empty() && tInfo.fUnsplittableMembers.empty()) continue;
         clang::ASTContext &C = recordDecl.getASTContext();
         clang::SourceRange commentRange; // Empty: this is a fake comment
         std::string userDefinedProperty;
         userDefinedProperty.reserve(100);
         for (auto fieldPtr : recordDecl.fields()) {
            const auto fieldName = fieldPtr->getName().data();
            if (tInfo.fTransientMembers.count(fieldName) == 1) {
               userDefinedProperty = "!";
            } else if (tInfo.fUnsplittableMembers.count(fieldName) == 1) {
               userDefinedProperty = propNames::comment + propNames::separator + "||";
            }
            if (!userDefinedProperty.empty()) {
               fieldPtr->addAttr(new(C) clang::AnnotateAttr(commentRange, C, userDefinedProperty, 0));
               userDefinedProperty = "";
            }
         }
      } // End loop on template info
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////

bool DictSelectionReader::VisitRecordDecl(clang::RecordDecl *recordDecl)
{
   if (auto CXXRD = llvm::dyn_cast<clang::CXXRecordDecl>(recordDecl)) {
      if (CXXRD->getDescribedClassTemplate()) {
         // this is a member of a template; ignore: dictionaries can only
         // select instances / specializations.
         return true;
      }
   }
   if (auto CXXRD = llvm::dyn_cast<clang::CXXRecordDecl>(recordDecl->getDeclContext())) {
      if (CXXRD->getDescribedClassTemplate()) {
         // this is a member of a template; ignore: dictionaries can only
         // select instances / specializations.
         return true;
      }
   }
   if (fIsFirstPass)
      return FirstPass(*recordDecl);
   else
      return SecondPass(*recordDecl);
}


////////////////////////////////////////////////////////////////////////////////

/**
 * Transform a name of a class instance into a pattern for selection
 * e.g. myClass<double, int, ...> in the selection namespace
 * will translate into a pattern of the type myClass<*>
 **/
inline std::string DictSelectionReader::PatternifyName(const std::string &className)
{
   return className.substr(0, className.find_first_of("<")) + "<*>";

}

////////////////////////////////////////////////////////////////////////////////

/**
 * Transform the name of the type eliminating the trailing & and *
 **/
inline void DictSelectionReader::GetPointeeType(std::string &typeName)
{
   while (typeName[typeName.size() - 1] == '*' ||
          typeName[typeName.size() - 1] == '&') {
      typeName = typeName.substr(0, typeName.size() - 1);
   }
}

}
}
