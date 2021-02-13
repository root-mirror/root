// @(#)root/core/utils:$Id: SelectionRules.cxx 41697 2011-11-01 21:03:41Z pcanal $
// Author: Velislava Spasova September 2010

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class SelectionRules
The class representing the collection of selection rules.
*/

#include <iostream>
#include <sstream>
#include <algorithm>
#ifndef WIN32
#include <fnmatch.h>
#else
#include "Shlwapi.h"
#define fnmatch(glob, path, dummy) PathMatchSpecA(path, glob);
#endif
#include "RtypesCore.h"
#include "SelectionRules.h"
#include "llvm/Support/raw_ostream.h"
#include "clang/Basic/SourceLocation.h"
#include "clang/Basic/SourceManager.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/DeclTemplate.h"

#include "cling/Interpreter/Interpreter.h"

const clang::CXXRecordDecl *R__ScopeSearch(const char *name, const clang::Type** resultType = 0) ;

void SelectionRules::AddClassSelectionRule(const ClassSelectionRule& classSel)
{
   fRulesCounter++;
   fClassSelectionRules.push_front(classSel);
   if (!classSel.HasInterpreter())
      fClassSelectionRules.begin()->SetInterpreter(fInterp);
   if (classSel.GetIndex() < 0)
      fClassSelectionRules.begin()->SetIndex(fRulesCounter);
}

void SelectionRules::AddFunctionSelectionRule(const FunctionSelectionRule& funcSel)
{
   fRulesCounter++;
   fFunctionSelectionRules.push_back(funcSel);
   if (!funcSel.HasInterpreter())
      fFunctionSelectionRules.begin()->SetInterpreter(fInterp);
   if (funcSel.GetIndex() < 0)
      fFunctionSelectionRules.begin()->SetIndex(fRulesCounter);
}

void SelectionRules::AddVariableSelectionRule(const  VariableSelectionRule& varSel)
{
   fRulesCounter++;
   fVariableSelectionRules.push_back(varSel);
   if (!varSel.HasInterpreter())
      fVariableSelectionRules.begin()->SetInterpreter(fInterp);
   if (varSel.GetIndex() < 0)
      fVariableSelectionRules.begin()->SetIndex(fRulesCounter);
}

void SelectionRules::AddEnumSelectionRule(const EnumSelectionRule& enumSel)
{
   fRulesCounter++;
   fEnumSelectionRules.push_back(enumSel);
   if (!enumSel.HasInterpreter())
      fEnumSelectionRules.begin()->SetInterpreter(fInterp);
   if (enumSel.GetIndex() < 0)
      fEnumSelectionRules.begin()->SetIndex( fRulesCounter );
}

void SelectionRules::PrintSelectionRules() const
{
   std::cout<<"Printing Selection Rules:"<<std::endl;
   if (!fClassSelectionRules.empty()) {
      int i = 0;
      for(std::list<ClassSelectionRule>::const_iterator it = fClassSelectionRules.begin();
          it != fClassSelectionRules.end(); ++it, ++i) {
         std::cout<<"\tClass sel rule "<<i<<":"<<std::endl;
         std::cout<< *it;
      }
   }
   else {
      std::cout<<"\tNo Class Selection Rules"<<std::endl;
   }

   if (!fFunctionSelectionRules.empty()) {
      //std::cout<<""<<std::endl;
      std::list<FunctionSelectionRule>::const_iterator it2;
      int i = 0;

      for (it2 = fFunctionSelectionRules.begin(); it2 != fFunctionSelectionRules.end(); ++it2, ++i) {
         std::cout<<"\tFunction sel rule "<<i<<":"<<std::endl;
         std::cout<<"\t\tSelected: ";
         switch(it2->GetSelected()){
            case BaseSelectionRule::kYes: std::cout<<"Yes"<<std::endl;
               break;
            case BaseSelectionRule::kNo: std::cout<<"No"<<std::endl;
               break;
            case BaseSelectionRule::kDontCare: std::cout<<"Don't Care"<<std::endl;
               break;
            default: std::cout<<"Unspecified"<<std::endl;
         }
         it2->PrintAttributes(std::cout,2);
      }
   }
   else {
      std::cout<<"\tNo function sel rules"<<std::endl;
   }

   if (!fVariableSelectionRules.empty()) {
      std::list<VariableSelectionRule>::const_iterator it3;
      int i = 0;

      for (it3 = fVariableSelectionRules.begin(); it3 != fVariableSelectionRules.end(); ++it3, ++i) {
         std::cout<<"\tVariable sel rule "<<i<<":"<<std::endl;
         std::cout<<"\t\tSelected: ";
         switch(it3->GetSelected()){
            case BaseSelectionRule::kYes: std::cout<<"Yes"<<std::endl;
               break;
            case BaseSelectionRule::kNo: std::cout<<"No"<<std::endl;
               break;
            case BaseSelectionRule::kDontCare: std::cout<<"Don't Care"<<std::endl;
               break;
            default: std::cout<<"Unspecified"<<std::endl;
         }
         it3->PrintAttributes(std::cout,2);
      }
   }
   else {
      std::cout<<"\tNo variable sel rules"<<std::endl;
   }

   if (!fEnumSelectionRules.empty()) {
      std::list<EnumSelectionRule>::const_iterator it4;
      int i = 0;

      for (it4 = fEnumSelectionRules.begin(); it4 != fEnumSelectionRules.end(); ++it4, ++i) {
         std::cout<<"\tEnum sel rule "<<i<<":"<<std::endl;
         std::cout<<"\t\tSelected: ";
         switch(it4->GetSelected()){
            case BaseSelectionRule::kYes: std::cout<<"Yes"<<std::endl;
               break;
            case BaseSelectionRule::kNo: std::cout<<"No"<<std::endl;
               break;
            case BaseSelectionRule::kDontCare: std::cout<<"Don't Care"<<std::endl;
               break;
            default: std::cout<<"Unspecified"<<std::endl;
         }
         it4->PrintAttributes(std::cout,2);
      }
   }
   else {
      std::cout<<"\tNo enum sel rules"<<std::endl;
   }
}

void SelectionRules::ClearSelectionRules()
{
   fClassSelectionRules.clear();
   fFunctionSelectionRules.clear();
   fVariableSelectionRules.clear();
   fEnumSelectionRules.clear();
}

template<class RULE>
static bool HasDuplicate(RULE* rule,
                         std::unordered_map<std::string,RULE*>& storedRules,
                         const std::string& attrName){
   auto itRetCodePair = storedRules.emplace( attrName, rule );

   auto storedRule = storedRules[attrName];

   if (itRetCodePair.second ||
       storedRule->GetSelected() != rule->GetSelected())  return false;
   auto areEqual = SelectionRulesUtils::areEqual(storedRule,rule);

   std::stringstream sstr; sstr << "Rule:\n";
   rule->Print(sstr);
   sstr << (areEqual ? "Identical " : "Conflicting ");
   sstr << "rule already stored:\n";
   storedRule->Print(sstr);
   ROOT::TMetaUtils::Warning("SelectionRules::CheckDuplicates",
                             "Duplicated rule found.\n%s",sstr.str().c_str());
   return !areEqual;
}

template<class RULESCOLLECTION, class RULE = typename RULESCOLLECTION::value_type>
static int CheckDuplicatesImp(RULESCOLLECTION& rules){
   int nDuplicates = 0;
   std::unordered_map<std::string, RULE*> patterns,names;
   for (auto&& rule : rules){
      if (rule.HasAttributeName() && HasDuplicate(&rule,names,rule.GetAttributeName())) nDuplicates++;
      if (rule.HasAttributePattern() && HasDuplicate(&rule,patterns,rule.GetAttributePattern())) nDuplicates++;
   }
   return nDuplicates;
}

int SelectionRules::CheckDuplicates(){

   int nDuplicates = 0;
   nDuplicates += CheckDuplicatesImp(fClassSelectionRules);
   nDuplicates += CheckDuplicatesImp(fFunctionSelectionRules);
   nDuplicates += CheckDuplicatesImp(fVariableSelectionRules);
   nDuplicates += CheckDuplicatesImp(fEnumSelectionRules);
   if (0 != nDuplicates){
      ROOT::TMetaUtils::Error("SelectionRules::CheckDuplicates",
            "Duplicates in rules were found.\n");
   }
   return nDuplicates;
}

static bool Implies(const ClassSelectionRule& patternRule, const ClassSelectionRule& nameRule){

   // Check if these both select or both exclude
   if (patternRule.GetSelected() != nameRule.GetSelected()) return false;

   // If the two rules are not compatible modulo their name/pattern, bail out
   auto nAttrsPattern = patternRule.GetAttributes().size();
   auto nAttrsName = nameRule.GetAttributes().size();
   if ((nAttrsPattern != 1 || nAttrsName !=1) &&
       !SelectionRulesUtils::areEqual(&patternRule, &nameRule, true /*moduloNameOrPattern*/)) {
      return false;
   }

   auto pattern = patternRule.GetAttributePattern().c_str();
   auto name = nameRule.GetAttributeName().c_str();

   // Now check if the pattern matches the name
   auto implies = 0 ==  fnmatch(pattern, name, FNM_PATHNAME);

   if (implies){
      static const auto msg = "The pattern rule %s matches the name rule %s. "
      "Since the name rule has compatible attributes, "
      "it will be removed: the pattern rule will match the necessary classes if needed.\n";

      ROOT::TMetaUtils::Info("SelectionRules::Optimize", msg, pattern, name);
   }


   return implies;

}

void SelectionRules::Optimize(){

   // Remove name rules "implied" by pattern rules

   if (!IsSelectionXMLFile()) return;

   const auto& selectionRules = fClassSelectionRules;

   auto predicate = [&selectionRules](const ClassSelectionRule &rule) -> bool {
     if (rule.HasAttributeName()) {
        for (auto&& intRule : selectionRules){
           if (intRule.HasAttributePattern() && Implies(intRule, rule)) {
              return true;
           }
        }
     }
     return false;
   };
   fClassSelectionRules.remove_if(predicate);
}

const ClassSelectionRule *SelectionRules::IsDeclSelected(const clang::RecordDecl *D, bool includeTypedefRule) const
{
   std::string qual_name;
   GetDeclQualName(D,qual_name);
   return IsClassSelected(D, qual_name, includeTypedefRule);
}

const ClassSelectionRule *SelectionRules::IsDeclSelected(const clang::TypedefNameDecl *D) const
{
   std::string qual_name;
   GetDeclQualName(D,qual_name);
   return IsClassSelected(D, qual_name, true);
}

const ClassSelectionRule *SelectionRules::IsDeclSelected(const clang::NamespaceDecl *D) const
{
   std::string qual_name;
   GetDeclQualName(D,qual_name);
   return IsNamespaceSelected(D, qual_name);
}

const BaseSelectionRule *SelectionRules::IsDeclSelected(const clang::EnumDecl *D) const
{
   // Currently rootcling does not need any information on enums, except
   // for the PCM / proto classes that register them to build TEnums without
   // parsing. This can be removed once (real) PCMs are available.
   // Note that the code below was *not* properly matching the case
   //   typedef enum { ... } abc;
   // as the typedef is stored as an anonymous EnumDecl in clang.
   // It is likely that using a direct lookup on the name would
   // return the appropriate typedef (and then we would need to
   // select 'both' the typedef and the anonymous enums.

#if defined(R__MUST_REVISIT)
# if R__MUST_REVISIT(6,4)
   "Can become no-op once PCMs are available."
# endif
#endif

   std::string str_name;   // name of the Decl
   std::string qual_name;  // fully qualified name of the Decl
   GetDeclName(D, str_name, qual_name);

   if (IsParentClass(D)) {
      const BaseSelectionRule *selector = IsMemberSelected(D, str_name);
      if (!selector) // if the parent class is deselected, we could still select the enum
         return IsEnumSelected(D, qual_name);
      else           // if the parent class is selected so are all nested enums
         return selector;
   }

   // Enum is not part of a class
   else {
      if (IsLinkdefFile())
         return IsLinkdefEnumSelected(D, qual_name);
      return IsEnumSelected(D, qual_name);
   }

   return 0;
}

const BaseSelectionRule *SelectionRules::IsDeclSelected(const clang::VarDecl* D) const
{
   std::string qual_name;  // fully qualified name of the Decl
   GetDeclQualName(D, qual_name);

   if (IsLinkdefFile())
      return IsLinkdefVarSelected(D, qual_name);
   else
      return IsVarSelected(D, qual_name);

}

const BaseSelectionRule *SelectionRules::IsDeclSelected(const clang::FieldDecl* /* D */) const
{
   // Currently rootcling does not need any information about fields.
   return 0;
#if 0
   std::string str_name;   // name of the Decl
   std::string qual_name;  // fully qualified name of the Decl
   GetDeclName(D, str_name, qual_name);

   return IsMemberSelected(D, str_name);
#endif
}

const BaseSelectionRule *SelectionRules::IsDeclSelected(const clang::FunctionDecl* D) const
{
   // Implement a simple matching for functions
   std::string qual_name;  // fully qualified name of the Decl
   GetDeclQualName(D, qual_name);
   if (IsLinkdefFile())
      return IsLinkdefFunSelected(D, qual_name);
   else
      return IsFunSelected(D, qual_name);
}

const BaseSelectionRule *SelectionRules::IsDeclSelected(const clang::Decl *D) const
{
   if (!D) {
      return 0;
   }

   clang::Decl::Kind declkind = D->getKind();

   switch (declkind) {
   case clang::Decl::CXXRecord:
   case clang::Decl::ClassTemplateSpecialization:
   case clang::Decl::ClassTemplatePartialSpecialization:
      // structs, unions and classes are all CXXRecords
      return IsDeclSelected(llvm::dyn_cast<clang::RecordDecl>(D));
   case clang::Decl::Namespace:
      return IsDeclSelected(llvm::dyn_cast<clang::NamespaceDecl>(D));
   case clang::Decl::Enum:
      // Enum is part of a class
      return IsDeclSelected(llvm::dyn_cast<clang::EnumDecl>(D));
   case clang::Decl::Var:
      return IsDeclSelected(llvm::dyn_cast<clang::VarDecl>(D));
#if ROOTCLING_NEEDS_FUNCTIONS_SELECTION
   case clang::Decl::Function:
      return IsDeclSelected(llvm::dyn_cast<clang::FunctionDecl>(D));
   case clang::Decl::CXXMethod:
   case clang::Decl::CXXConstructor:
   case clang::Decl::CXXDestructor: {
      // std::string proto;
      //  if (GetFunctionPrototype(D,proto))
      //       std::cout<<"\n\tFunction prototype: "<<str_name + proto;
      //  else
      //       std::cout<<"Error in prototype formation";
      std::string str_name;   // name of the Decl
      std::string qual_name;  // fully qualified name of the Decl
      GetDeclName(D, str_name, qual_name);
      if (IsLinkdefFile()) {
         return IsLinkdefMethodSelected(D, qual_name);
      }
      return IsMemberSelected(D, str_name);
   }
#endif
   case clang::Decl::Field:
      return IsDeclSelected(llvm::dyn_cast<clang::FieldDecl>(D));
   default:
      // Humm we are not treating this case!
      return 0;
   }

   // std::string str_name;   // name of the Decl
   // std::string qual_name;  // fully qualified name of the Decl
   // GetDeclName(D, str_name, qual_name);
   // fprintf(stderr,"IsDeclSelected: %s %s\n", str_name.c_str(), qual_name.c_str());
}


bool SelectionRules::GetDeclName(const clang::Decl* D, std::string& name, std::string& qual_name) const
{
   const clang::NamedDecl* N = llvm::dyn_cast<clang::NamedDecl> (D);

   if (!N)
      return false;

   // the identifier is NULL for some special methods like constructors, destructors and operators
   if (N->getIdentifier()) {
      name = N->getNameAsString();
   }
   else if (N->isCXXClassMember()) { // for constructors, destructors, operator=, etc. methods
      name =  N->getNameAsString(); // we use this (unefficient) method to Get the name in that case
   }
   llvm::raw_string_ostream stream(qual_name);
   N->getNameForDiagnostic(stream,N->getASTContext().getPrintingPolicy(),true);
   return true;
}

void SelectionRules::GetDeclQualName(const clang::Decl* D, std::string& qual_name) const{
   const clang::NamedDecl* N = static_cast<const clang::NamedDecl*> (D);
   llvm::raw_string_ostream stream(qual_name);
   N->getNameForDiagnostic(stream,N->getASTContext().getPrintingPolicy(),true);
   }

bool SelectionRules::GetFunctionPrototype(const clang::FunctionDecl* F, std::string& prototype) const {

   if (!F) {
      return false;
   }

   const std::vector<std::string> quals={"*","&"};

   prototype = "";
   // iterate through all the function parameters
   std::string type;
   for (auto I = F->param_begin(), E = F->param_end(); I != E; ++I) {

      clang::ParmVarDecl* P = *I;

      if (prototype != "")
         prototype += ",";

      ROOT::TMetaUtils::GetNormalizedName(type,P->getType(),fInterp,fNormCtxt);

      // We need to get rid of the "class " string if present
      ROOT::TMetaUtils::ReplaceAll(type,"class ", "");
      // We need to get rid of the "restrict " string if present
      ROOT::TMetaUtils::ReplaceAll(type,"restrict", "");

      // pointers are returned in the form "int *" and I need them in the form "int*"
      // same for &
      for (auto& qual : quals){
        auto pos = type.find(" "+qual);
        if (pos != std::string::npos)
           type.replace( pos, 2, qual );
        }
//          for (auto& qual : quals){
//             if (type.at(type.length()-1) == qual && type.at(type.length()-2) == ' ') {
//                type.at(type.length()-2) = qual;
//                type.erase(type.length()-1);
//             }
//          }
      prototype += type;
   }
   prototype = "(" + prototype + ")";
   return true;

}


bool SelectionRules::IsParentClass(const clang::Decl* D) const
{
   //TagDecl has methods to understand of what kind is the Decl - class, struct or union
   if (const clang::TagDecl* T
       = llvm::dyn_cast<clang::TagDecl>(D->getDeclContext()))
      return T->isClass() || T->isStruct();
   return false;
}


bool SelectionRules::IsParentClass(const clang::Decl* D, std::string& parent_name, std::string& parent_qual_name) const
{
   if (const clang::TagDecl* parent
       = llvm::dyn_cast<clang::TagDecl>(D->getDeclContext())) {
      if (parent->isClass()|| parent->isStruct()) {
         GetDeclName(parent, parent_name, parent_qual_name);
         return true;
      }
   }
   return false;
}

bool SelectionRules::GetParentName(const clang::Decl* D, std::string& parent_name, std::string& parent_qual_name) const
{
   if (const clang::RecordDecl* parent
       = llvm::dyn_cast<clang::RecordDecl>(D->getDeclContext())) {
      GetDeclName(parent, parent_name, parent_qual_name);
      return true;
   }
   return false;
}

/* This is the method that crashes
 bool SelectionRules::GetParent(clang::Decl* D, clang::Decl* parent)
 {
 clang::DeclContext *ctx = D->GetDeclContext();

 if (ctx->isRecord()){
 //DEBUG std::cout<<"\n\tDeclContext is Record";
 parent = llvm::dyn_cast<clang::Decl> (ctx);
 if (!parent) {
 return false;
 }
 else {
 return true;
 }
 }
 else return false;
 }
 */


// isClassSelected checks if a class is selected or not. Thre is a difference between the
// behaviour of rootcint and genreflex especially with regard to class pattern processing.
// In genreflex if we have <class pattern = "*" /> this will select all the classes
// (and structs) found in the header file. In rootcint if we have something similar, i.e.
// #pragma link C++ class *, we will select only the outer classes - for the nested
// classes we have to specifie #pragma link C++ class *::*. And yet this is only valid
// for one level of nesting - if you need it for many levels of nesting, you will
// probably have to implement it yourself.
// Here the idea is the following - we traverse the list of class selection rules.
// For every class we check do we have a class selection rule. We use here the
// method isSelected() (defined in BaseSelectionRule.cxx). This method returns true
// only if we have class selection rule which says "Select". Otherwise it returns
// false - in which case we have to check wether we found a class selection rule
// which says "Veto" (noName = false and don't Care = false; OR noName = false and
// don't Care = true and we don't have neither method nor field selection rules -
// which is for the selection.xml file case). If noName is true than we just continue -
// this means that the current class selection rule isn't applicable for this class.

const ClassSelectionRule *SelectionRules::IsNamespaceSelected(const clang::Decl* D, const std::string& qual_name) const
{
   const clang::NamespaceDecl* N = llvm::dyn_cast<clang::NamespaceDecl> (D); //TagDecl has methods to understand of what kind is the Decl
   if (N==0) {
      std::cout<<"\n\tCouldn't cast Decl to NamespaceDecl";
      return 0;
   }

   const ClassSelectionRule *selector = 0;
   int fImplNo = 0;
   const ClassSelectionRule *explicit_selector = 0;
   const ClassSelectionRule *specific_pattern_selector = 0;
   int fFileNo = 0;

   // NOTE: should we separate namespaces from classes in the rules?
   std::list<ClassSelectionRule>::const_iterator it = fClassSelectionRules.begin();
   // iterate through all class selection rles
   std::string name_value;
   std::string pattern_value;
   BaseSelectionRule::EMatchType match;
   for(; it != fClassSelectionRules.end(); ++it) {

      match = it->Match(N,qual_name,"",IsLinkdefFile());

      if (match != BaseSelectionRule::kNoMatch) {
         // If we have a match.
         if (it->GetSelected() == BaseSelectionRule::kYes) {
            selector = &(*it);
            if (IsLinkdefFile()){
               // rootcint prefers explicit rules over pattern rules
               if (match == BaseSelectionRule::kName) {
                  explicit_selector = &(*it);
               } else if (match == BaseSelectionRule::kPattern) {
                  // NOTE: weird ...
                  if (it->GetAttributeValue("pattern", pattern_value) &&
                      pattern_value != "*" && pattern_value != "*::*") specific_pattern_selector = &(*it);
               }
            }
         } else if (it->GetSelected() == BaseSelectionRule::kNo) {
            if (!IsLinkdefFile()) {
               // in genreflex - we could explicitly select classes from other source files
               if (match == BaseSelectionRule::kFile) ++fFileNo; // if we have veto because of class defined in other source file -> implicit No
               else {

#ifdef SELECTION_DEBUG
                  std::cout<<"\tNo returned"<<std::endl;
#endif

                  return 0; // explicit No returned
               }
            }
            if (match == BaseSelectionRule::kPattern) {
               //this is for the Linkdef selection
               if (it->GetAttributeValue("pattern", pattern_value) &&
                   (pattern_value == "*" || pattern_value == "*::*")) ++fImplNo;
               else
                  return 0;
            }
            else
               return 0;
         }
         else if (it->GetSelected() == BaseSelectionRule::kDontCare && !(it->HasMethodSelectionRules()) && !(it->HasFieldSelectionRules())) {

#ifdef SELECTION_DEBUG
            std::cout<<"Empty dontC returned = No"<<std::endl;
#endif

            return 0;
         }
      }
   }
   if (IsLinkdefFile()) {
      // for rootcint explicit (name) Yes is stronger than implicit (pattern) No which is stronger than implicit (pattern) Yes

#ifdef SELECTION_DEBUG
      std::cout<<"\n\tfYes = "<<fYes<<", fImplNo = "<<fImplNo<<std::endl;
#endif

      if (explicit_selector) return explicit_selector;
      else if (specific_pattern_selector) return specific_pattern_selector;
      else if (fImplNo > 0) return 0;
      else return selector;
   }
   else {
      // for genreflex explicit Yes is stronger than implicit file No

#ifdef SELECTION_DEBUG
      std::cout<<"\n\tfYes = "<<fYes<<", fFileNo = "<<fFileNo<<std::endl;
#endif

      if (selector)
         return selector;
      else
         return 0;
   }

}


const ClassSelectionRule *SelectionRules::IsClassSelected(const clang::Decl* D, const std::string& qual_name, bool includeTypedefRule) const
{
   const clang::TagDecl* tagDecl = llvm::dyn_cast<clang::TagDecl> (D); //TagDecl has methods to understand of what kind is the Decl
   const clang::TypedefNameDecl* typeDefNameDecl = llvm::dyn_cast<clang::TypedefNameDecl> (D);

   if (!tagDecl && !typeDefNameDecl) { // Ill posed
      ROOT::TMetaUtils::Error("SelectionRules::IsClassSelected",
            "Cannot cast Decl to TagDecl and Decl is not a typedef.\n");
      return 0;
      }

   if (!tagDecl && typeDefNameDecl){ // Let's try to fetch the underlying RecordDecl
      clang::RecordDecl* recordDecl = ROOT::TMetaUtils::GetUnderlyingRecordDecl(typeDefNameDecl->getUnderlyingType());
      if (!recordDecl){
         ROOT::TMetaUtils::Error("SelectionRules::IsClassSelected",
                                 "Cannot get RecordDecl behind TypedefDecl.\n");
         return 0;
      }
      tagDecl = recordDecl;
   }

   // At this point tagDecl must be well defined
   const bool isLinkDefFile =  IsLinkdefFile();
   if (!( isLinkDefFile || tagDecl->isClass() || tagDecl->isStruct() ))
      return 0; // Union for Genreflex

   const ClassSelectionRule *selector = 0;
   int fImplNo = 0;
   const ClassSelectionRule *explicit_selector = 0;
   const ClassSelectionRule *specific_pattern_selector = 0;
   int fFileNo = 0;

   // iterate through all class selection rles
   bool earlyReturn=false;
   const ClassSelectionRule* retval = nullptr;
   const clang::NamedDecl* nDecl(llvm::dyn_cast<clang::NamedDecl>(D));
   for(auto& rule : fClassSelectionRules) {
      if (!includeTypedefRule && rule.IsFromTypedef())
         continue;
      BaseSelectionRule::EMatchType match = rule.Match(nDecl, qual_name, "", isLinkDefFile);
      if (match != BaseSelectionRule::kNoMatch) {
         // Check if the template must have its arguments manipulated
         if (const clang::ClassTemplateSpecializationDecl* ctsd =
         llvm::dyn_cast_or_null<clang::ClassTemplateSpecializationDecl>(D))
            if(const clang::ClassTemplateDecl* ctd = ctsd->getSpecializedTemplate()){
               const std::string& nArgsToKeep = rule.GetAttributeNArgsToKeep();
               if (!nArgsToKeep.empty()){
                  fNormCtxt.AddTemplAndNargsToKeep(ctd->getCanonicalDecl(),
                                                   std::atoi(nArgsToKeep.c_str()));
               }
            }

         if (earlyReturn) continue;

         // If we have a match.
         selector = &(rule);
         if (rule.GetSelected() == BaseSelectionRule::kYes) {

            if (isLinkDefFile){
               // rootcint prefers explicit rules over pattern rules
               if (match == BaseSelectionRule::kName) {
                  explicit_selector = &(rule);
               } else if (match == BaseSelectionRule::kPattern) {
                  // NOTE: weird ...
                  const std::string& pattern_value=rule.GetAttributePattern();
                  if (!pattern_value.empty() &&
                      pattern_value != "*" &&
                      pattern_value != "*::*") specific_pattern_selector = &(rule);
               }
            }
         } else if (rule.GetSelected() == BaseSelectionRule::kNo) {

            if (!isLinkDefFile) {
               // in genreflex - we could explicitly select classes from other source files
               if (match == BaseSelectionRule::kFile) ++fFileNo; // if we have veto because of class defined in other source file -> implicit No
               else {
                  retval = selector;
                  earlyReturn=true; // explicit No returned
               }
            }
            if (match == BaseSelectionRule::kPattern) {
               //this is for the Linkdef selection
               const std::string& pattern_value=rule.GetAttributePattern();
               if (!pattern_value.empty() &&
                   (pattern_value == "*" || pattern_value == "*::*")) ++fImplNo;
               else
                  earlyReturn=true;
            }
            else
               earlyReturn=true;
         }
         else if (rule.GetSelected() == BaseSelectionRule::kDontCare && !(rule.HasMethodSelectionRules()) && !(rule.HasFieldSelectionRules())) {
            earlyReturn=true;
         }
      }
   } // Loop over the rules.

   if (earlyReturn) return retval;

   if (isLinkDefFile) {
      // for rootcint explicit (name) Yes is stronger than implicit (pattern) No which is stronger than implicit (pattern) Yes
      if (explicit_selector) return explicit_selector;
      else if (specific_pattern_selector) return specific_pattern_selector;
      else if (fImplNo > 0) return 0;
      else return selector;
   }
   else {
      // for genreflex explicit Yes is stronger than implicit file No
      return selector; // it can be nullptr
   }

}

const BaseSelectionRule *SelectionRules::IsVarSelected(const clang::VarDecl* D, const std::string& qual_name) const
{
   std::list<VariableSelectionRule>::const_iterator it = fVariableSelectionRules.begin();
   std::list<VariableSelectionRule>::const_iterator it_end =  fVariableSelectionRules.end();

   const BaseSelectionRule *selector = 0;

   // iterate through all the rules
   // we call this method only for genrefex variables, functions and enums - it is simpler than the class case:
   // if we have No - it is veto even if we have explicit yes as well
   for(; it != it_end; ++it) {
      if (BaseSelectionRule::kNoMatch != it->Match(llvm::dyn_cast<clang::NamedDecl>(D), qual_name, "", false)) {
         if (it->GetSelected() == BaseSelectionRule::kNo) {
            // The rule did explicitly request to not select this entity.
            return 0;
         } else {
            selector = &(*it);
         }
      }
   }

   return selector;
}

const BaseSelectionRule *SelectionRules::IsFunSelected(const clang::FunctionDecl *D, const std::string &qual_name) const
{

   if (fFunctionSelectionRules.size() == 0 ||
       D->getPrimaryTemplate() != nullptr ||
       llvm::isa<clang::CXXMethodDecl>(D)) return nullptr;

   std::string prototype;
   GetFunctionPrototype(D, prototype);
   prototype = qual_name + prototype;

   const BaseSelectionRule *selector = nullptr;
   // iterate through all the rules
   // we call this method only for genrefex variables, functions and enums - it is simpler than the class case:
   // if we have No - it is veto even if we have explicit yes as well
   for (const auto & rule : fFunctionSelectionRules) {
      if (BaseSelectionRule::kNoMatch != rule.Match(D, qual_name, prototype, false)) {
         if (rule.GetSelected() == BaseSelectionRule::kNo) {
            // The rule did explicitly request to not select this entity.
            return nullptr;
         } else {
            selector = &(rule);
         }
      }
   }

   return selector;
}

const BaseSelectionRule *SelectionRules::IsEnumSelected(const clang::EnumDecl* D, const std::string& qual_name) const
{
   const BaseSelectionRule *selector = 0;

   // iterate through all the rules
   // we call this method only for genrefex variables, functions and enums - it is simpler than the class case:
   // if we have No - it is veto even if we have explicit yes as well
   for(const auto& rule: fEnumSelectionRules) {
      if (BaseSelectionRule::kNoMatch != rule.Match(llvm::dyn_cast<clang::NamedDecl>(D), qual_name, "", false)) {
         if (rule.GetSelected() == BaseSelectionRule::kNo) {
            // The rule did explicitly request to not select this entity.
            return 0;
         } else {
            selector = &rule;
         }
      }
   }

   return selector;
}

const BaseSelectionRule *SelectionRules::IsLinkdefVarSelected(const clang::VarDecl* D, const std::string& qual_name) const
{


   const BaseSelectionRule *selector = 0;
   int fImplNo = 0;
   const BaseSelectionRule *explicit_selector = 0;

   std::string name_value;
   std::string pattern_value;
   for(auto& selRule: fVariableSelectionRules) {
      BaseSelectionRule::EMatchType match
         = selRule.Match(llvm::dyn_cast<clang::NamedDecl>(D), qual_name, "", false);
      if (match != BaseSelectionRule::kNoMatch) {
         if (selRule.GetSelected() == BaseSelectionRule::kYes) {
            // explicit rules are with stronger priority in rootcint
            if (IsLinkdefFile()){
               if (match == BaseSelectionRule::kName) {
                  explicit_selector = &selRule;
               } else if (match == BaseSelectionRule::kPattern) {
                  if (selRule.GetAttributeValue("pattern", pattern_value)) {
                     explicit_selector=&selRule;
                     // NOTE: Weird ... This is a strange definition of explicit.
                     //if (pattern_value != "*" && pattern_value != "*::*") explicit_selector = selRule;
                  }
               }
            }
         }
         else {
            if (!IsLinkdefFile()) return 0;
            else {
               if (selRule.GetAttributeValue("pattern", pattern_value)) {
                  if (pattern_value == "*" || pattern_value == "*::*") ++fImplNo;
                  else
                     return 0;
               }
               else
                  return 0;
            }
         }
      }
   }

   if (IsLinkdefFile()) {

#ifdef SELECTION_DEBUG
      std::cout<<"\n\tfYes = "<<fYes<<", fImplNo = "<<fImplNo<<std::endl;
#endif

      if (explicit_selector) return explicit_selector;
      else if (fImplNo > 0) return 0;
      else return selector;
   }
   else{
      return selector;
   }
}

const BaseSelectionRule *SelectionRules::IsLinkdefFunSelected(const clang::FunctionDecl* D, const std::string& qual_name) const
{

   if (fFunctionSelectionRules.size() == 0 ||
       D->getPrimaryTemplate() != nullptr ||
       llvm::isa<clang::CXXMethodDecl>(D)) return nullptr;

   std::string prototype;

   GetFunctionPrototype(D, prototype);
   prototype = qual_name + prototype;

   const BaseSelectionRule *selector = 0;
   int fImplNo = 0;
   const BaseSelectionRule *explicit_selector = 0;

   std::string pattern_value;
   for(auto& selRule : fFunctionSelectionRules) {
      BaseSelectionRule::EMatchType match
         = selRule.Match(llvm::dyn_cast<clang::NamedDecl>(D), qual_name, prototype, false);
      if (match != BaseSelectionRule::kNoMatch) {
         if (selRule.GetSelected() == BaseSelectionRule::kYes) {
            // explicit rules are with stronger priority in rootcint
            if (IsLinkdefFile()){
               if (match == BaseSelectionRule::kName) {
                  explicit_selector = &selRule;
               } else if (match == BaseSelectionRule::kPattern) {
                  if (selRule.GetAttributeValue("pattern", pattern_value)) {
                     explicit_selector = &selRule;
                     // NOTE: Weird ... This is a strange definition of explicit.
                     //if (pattern_value != "*" && pattern_value != "*::*") explicit_selector = &selRule;
                  }
               }
            }
         }
         else {
            if (!IsLinkdefFile()) return 0;
            else {
               if (selRule.GetAttributeValue("pattern", pattern_value)) {
                  if (pattern_value == "*" || pattern_value == "*::*") ++fImplNo;
                  else
                     return 0;
               }
               else
                  return 0;
            }
         }
      }
   }

   if (IsLinkdefFile()) {
      if (explicit_selector) return explicit_selector;
      else if (fImplNo > 0) return 0;
      else return selector;
   }
   else{
      return selector;
   }
}

const BaseSelectionRule *SelectionRules::IsLinkdefEnumSelected(const clang::EnumDecl* D, const std::string& qual_name) const
{
   std::list<VariableSelectionRule>::const_iterator it;
   std::list<VariableSelectionRule>::const_iterator it_end;

   it = fEnumSelectionRules.begin();
   it_end = fEnumSelectionRules.end();

   const BaseSelectionRule *selector = 0;
   int fImplNo = 0;
   const BaseSelectionRule *explicit_selector = 0;

   std::string name_value;
   std::string pattern_value;
   for(; it != it_end; ++it) {
      BaseSelectionRule::EMatchType match =
         it->Match(llvm::dyn_cast<clang::NamedDecl>(D), qual_name, "", false);
      if (match != BaseSelectionRule::kNoMatch) {
         if (it->GetSelected() == BaseSelectionRule::kYes) {
            // explicit rules are with stronger priority in rootcint
            if (IsLinkdefFile()){
               if (match == BaseSelectionRule::kName){
                  explicit_selector = &(*it);
               } else if (match == BaseSelectionRule::kPattern &&
                          it->GetAttributeValue("pattern", pattern_value)) {
                  // Note: Weird ... This is a strange definition of explicit.
                  if (pattern_value != "*" && pattern_value != "*::*") explicit_selector = &(*it);
               }
            }
         }
         else {
            if (!IsLinkdefFile()) return 0;
            else {
               if (it->GetAttributeValue("pattern", pattern_value)) {
                  if (pattern_value == "*" || pattern_value == "*::*") ++fImplNo;
                  else
                     return 0;
               }
               else
                  return 0;
            }
         }
      }
   }

   if (IsLinkdefFile()) {

#ifdef SELECTION_DEBUG
      std::cout<<"\n\tfYes = "<<fYes<<", fImplNo = "<<fImplNo<<std::endl;
#endif

      if (explicit_selector) return explicit_selector;
      else if (fImplNo > 0) return 0;
      else return selector;
   }
   else{
      return selector;
   }
}

// In rootcint we could select and deselect methods independantly of the class/struct/union rules
// That's why we first have to check the explicit rules for the functions - to see if there
// is rule corresponding to our method.
// Which is more - if we have (and we can have) a pattern for the parent class, than a pattern for the
// nested class, than a pattern for certain methods in the nested class, than a rule for a
// method (name or prototype) in the nested class - the valid rule is the last one.
// This is true irrespective of the rules (select/deselect). This is not the case for genreflex -
// in genreflex if there is a pattern deselecting something even if we have an explicit rule
// to select it, it still will not be selected.
// This method (isLinkdefMethodSelected()) might be incomplete (but I didn't have the time to think
// of anything better)
//

const BaseSelectionRule *SelectionRules::IsLinkdefMethodSelected(const clang::Decl* D, const std::string& qual_name) const
{
   std::list<FunctionSelectionRule>::const_iterator it = fFunctionSelectionRules.begin();
   std::list<FunctionSelectionRule>::const_iterator it_end = fFunctionSelectionRules.end();
   std::string prototype;

   if (const clang::FunctionDecl* F = llvm::dyn_cast<clang::FunctionDecl> (D))
      GetFunctionPrototype(F, prototype);
   prototype = qual_name + prototype;

#ifdef SELECTION_DEBUG
   std::cout<<"\tFunction prototype = "<<prototype<<std::endl;
#endif

   int expl_Yes = 0, impl_r_Yes = 0, impl_rr_Yes = 0;
   int impl_r_No = 0, impl_rr_No = 0;
   const BaseSelectionRule *explicit_r = 0;
   const BaseSelectionRule *implicit_r = 0;
   const BaseSelectionRule *implicit_rr = 0;

   if (D->getKind() == clang::Decl::CXXMethod){
      // we first check the explicit rules for the method (in case of constructors and destructors we check the parent)
      std::string pat_value;
      for(; it != it_end; ++it) {
         BaseSelectionRule::EMatchType match
            = it->Match(llvm::dyn_cast<clang::NamedDecl>(D), qual_name, prototype, false);

         if (match == BaseSelectionRule::kName) {
            // here I should implement my implicit/explicit thing
            // I have included two levels of implicitness - "A::Get_*" is stronger than "*"
            explicit_r = &(*it);
            if (it->GetSelected() == BaseSelectionRule::kYes) ++expl_Yes;
            else {

#ifdef SELECTION_DEBUG
                  std::cout<<"\tExplicit rule BaseSelectionRule::kNo found"<<std::endl;
#endif

                  return 0; // == explicit BaseSelectionRule::kNo

            }
         } else if (match == BaseSelectionRule::kPattern) {

            if (it->GetAttributeValue("pattern", pat_value)) {
               if (pat_value == "*") continue; // we discard the global selection rules

               std::string par_name, par_qual_name;
               GetParentName(D, par_name, par_qual_name);
               std::string par_pat = par_qual_name + "::*";

               if (pat_value == par_pat) {
                  implicit_rr = &(*it);
                  if (it->GetSelected() == BaseSelectionRule::kYes) {

#ifdef SELECTION_DEBUG
                     std::cout<<"Implicit_rr rule ("<<pat_value<<"), selected = "<<selected<<std::endl;
#endif

                     ++impl_rr_Yes;
                  }
                  else {

#ifdef SELECTION_DEBUG
                     std::cout<<"Implicit_rr rule ("<<pat_value<<"), selected = "<<selected<<std::endl;
#endif

                     ++impl_rr_No;
                  }
               }
               else {
                  implicit_r = &(*it);
                  if (it->GetSelected() == BaseSelectionRule::kYes) {

#ifdef SELECTION_DEBUG
                     std::cout<<"Implicit_r rule ("<<pat_value<<"), selected = "<<selected<<std::endl;
#endif

                     ++impl_r_Yes;
                  }
                  else {

#ifdef SELECTION_DEBUG
                     std::cout<<"Implicit_r rule ("<<pat_value<<"), selected = "<<selected<<std::endl;
#endif

                     ++impl_r_No;
                  }
               }
            }
         }
      }
   }
   if (explicit_r /*&& expl_Yes > 0*/){

#ifdef SELECTION_DEBUG
      std::cout<<"\tExplicit rule BaseSelectionRule::BaseSelectionRule::kYes found"<<std::endl;
#endif

      return explicit_r; // if we have excplicit BaseSelectionRule::kYes
   }
   else if (implicit_rr) {
      if (impl_rr_No > 0) {

#ifdef SELECTION_DEBUG
         std::cout<<"\tImplicit_rr rule BaseSelectionRule::kNo found"<<std::endl;
#endif

         return 0;
      }
      else {

#ifdef SELECTION_DEBUG
         std::cout<<"\tImplicit_rr rule BaseSelectionRule::kYes found"<<std::endl;
#endif

         return implicit_rr;
      }
   }
   else if (implicit_r) {
      if (impl_r_No > 0) {

#ifdef SELECTION_DEBUG
         std::cout<<"\tImplicit_r rule BaseSelectionRule::kNo found"<<std::endl;
#endif

         return 0;
      }
      else {

#ifdef SELECTION_DEBUG
         std::cout<<"\tImplicit_r rule BaseSelectionRule::kYes found"<<std::endl;
#endif

         return implicit_r;
      }
   }
   else {

#ifdef SELECTION_DEBUG
      std::cout<<"\tChecking parent class rules"<<std::endl;
#endif
      // check parent


      std::string parent_name, parent_qual_name;
      if (!GetParentName(D, parent_name, parent_qual_name)) return 0;

      const BaseSelectionRule *selector = 0;
      int fImplNo = 0;
      const BaseSelectionRule *explicit_selector = 0;

      // the same as with isClass selected
      // I wanted to use GetParentDecl and then to pass is to isClassSelected because I didn't wanted to repeat
      // code but than GetParentDecl crashes (or returns non-sence Decl) for the built-in structs (__va_*)
      std::list<ClassSelectionRule>::const_iterator it = fClassSelectionRules.begin();
      std::string name_value;
      std::string pattern_value;
      for(; it != fClassSelectionRules.end(); ++it) {
         BaseSelectionRule::EMatchType match
            = it->Match(llvm::dyn_cast<clang::NamedDecl>(D), parent_qual_name, "", true); // == BaseSelectionRule::kYes

         if (match != BaseSelectionRule::kNoMatch) {
            if (it->GetSelected() == BaseSelectionRule::kYes) {
               selector = &(*it);

               if (match == BaseSelectionRule::kName) {
                  explicit_selector = &(*it);
               } else if (match == BaseSelectionRule::kPattern) {
                  if (it->GetAttributeValue("pattern", pattern_value)) {
                     // NOTE: weird ...
                     if (pattern_value != "*" && pattern_value != "*::*") explicit_selector = &(*it);
                  }
               }
            }
            else { // == BaseSelectionRule::kNo

               if (it->GetAttributeValue("pattern", pattern_value)) {
                  if (pattern_value == "*" || pattern_value == "*::*") ++fImplNo;
                  else
                     return 0;
               }
               else
                  return 0;
            }
         }
      }

#ifdef SELECTION_DEBUG
      std::cout<<"\n\tfYes = "<<fYes<<", fImplNo = "<<fImplNo<<std::endl;
#endif

      if (explicit_selector) {

#ifdef SELECTION_DEBUG
         std::cout<<"\tReturning Yes"<<std::endl;
#endif

         return explicit_selector;
      }
      else if (fImplNo > 0) {
#ifdef SELECTION_DEBUG
         std::cout<<"\tReturning No"<<std::endl;
#endif

         return 0;
      }
      else {

#ifdef SELECTION_DEBUG
         std::cout<<"\tReturning Yes"<<std::endl;
#endif

         return selector;
      }
   }

   return 0;

}

const BaseSelectionRule *SelectionRules::IsMemberSelected(const clang::Decl* D, const std::string& str_name) const
{
   std::string parent_name;
   std::string parent_qual_name;

   if (IsParentClass(D))
   {
      if (!GetParentName(D, parent_name, parent_qual_name)) return 0;

      const BaseSelectionRule *selector = 0;
      Int_t fImplNo = 0;
      const BaseSelectionRule *explicit_selector = 0;
      int fFileNo = 0;

      //DEBUG std::cout<<"\n\tParent is class";
      std::list<ClassSelectionRule>::const_iterator it = fClassSelectionRules.begin();
      std::string name_value;
      std::string pattern_value;
      for(; it != fClassSelectionRules.end(); ++it) {
         BaseSelectionRule::EMatchType match
            = it->Match(llvm::dyn_cast<clang::NamedDecl>(D), parent_qual_name, "", false);

         if (match != BaseSelectionRule::kNoMatch) {
            if (it->GetSelected() == BaseSelectionRule::kYes) {
               selector = &(*it);
               if (IsLinkdefFile()) {
                  if (match == BaseSelectionRule::kName) {
                     explicit_selector = &(*it);
                  } else if (match == BaseSelectionRule::kPattern) {
                     if (it->GetAttributeValue("pattern", pattern_value)) {
                        // NOTE: Weird ... This is a strange definition of explicit.
                        if (pattern_value != "*" && pattern_value != "*::*") explicit_selector = &(*it);
                     }
                  }
               }
            } else if (it->GetSelected() == BaseSelectionRule::kNo) {
               if (!IsLinkdefFile()) {
                  if (match == BaseSelectionRule::kFile) ++fFileNo;
                  else {

#ifdef SELECTION_DEBUG
                     std::cout<<"\tNo returned"<<std::endl;
#endif

                     return 0; // in genreflex we can't have that situation
                  }
               }
               else {
                  if (match == BaseSelectionRule::kPattern && it->GetAttributeValue("pattern", pattern_value)) {
                     if (pattern_value == "*" || pattern_value == "*::*") ++fImplNo;
                     else
                        return 0;
                  }
                  else
                     return 0;
               }
            }
            else if (it->GetSelected() == BaseSelectionRule::kDontCare) // - we check the method and field selection rules for the class
            {
               if (!it->HasMethodSelectionRules() && !it->HasFieldSelectionRules()) {

#ifdef SELECTION_DEBUG
                  std::cout<<"\tNo fields and methods"<<std::endl;
#endif

                  return 0; // == BaseSelectionRule::kNo
               }
               else {
                  clang::Decl::Kind kind = D->getKind();
                  if (kind == clang::Decl::Field || kind == clang::Decl::CXXMethod || kind == clang::Decl::CXXConstructor || kind == clang::Decl::CXXDestructor){
                     std::list<VariableSelectionRule> members;
                     std::list<VariableSelectionRule>::iterator mem_it;
                     std::list<VariableSelectionRule>::iterator mem_it_end;
                     std::string prototype;

                     if (kind == clang::Decl::Field) {
                        members = it->GetFieldSelectionRules();
                     }
                     else {
                        if (const clang::FunctionDecl* F = llvm::dyn_cast<clang::FunctionDecl> (D)){
                           GetFunctionPrototype(F, prototype);
                           prototype = str_name + prototype;
                        }
                        else{
                           ROOT::TMetaUtils::Warning("","Warning: could not cast Decl to FunctionDecl\n");
                        }
                        members = it->GetMethodSelectionRules();
                     }
                     mem_it = members.begin();
                     mem_it_end = members.end();
                     for (; mem_it != mem_it_end; ++mem_it) {
                        if (BaseSelectionRule::kName == mem_it->Match(llvm::dyn_cast<clang::NamedDecl>(D), str_name, prototype, false)) {
                           if (mem_it->GetSelected() == BaseSelectionRule::kNo) return 0;
                        }
                     }
                  }
               }
            }
         }
      }

      if (IsLinkdefFile()) {

#ifdef SELECTION_DEBUG
         std::cout<<"\n\tfYes = "<<fYes<<", fImplNo = "<<fImplNo<<std::endl;
#endif

         if (explicit_selector) {
#ifdef SELECTION_DEBUG
            std::cout<<"\tReturning Yes"<<std::endl;
#endif

            return explicit_selector;
         }
         else if (fImplNo > 0) {

#ifdef SELECTION_DEBUG
            std::cout<<"\tReturning No"<<std::endl;
#endif

            return 0;
         }
         else {

#ifdef SELECTION_DEBUG
            std::cout<<"\tReturning Yes"<<std::endl;
#endif

            return selector;
         }
      }
      else {

         return selector;
      }
   }
   else {
      return 0;
   }
}

bool SelectionRules::AreAllSelectionRulesUsed() const {
   for(auto&& rule : fClassSelectionRules){
      if (BaseSelectionRule::kNo!=rule.GetSelected() && !rule.GetMatchFound() /* && !GetHasFileNameRule() */ ) {
         std::string name;
         if (rule.GetAttributeValue("pattern", name)) {
            // keep it
         } else if (rule.GetAttributeValue("name", name)) {
            // keept it
         } else {
            name.clear();
         }
         std::string file_name_value;
         if (!rule.GetAttributeValue("file_name", file_name_value)) file_name_value.clear();

         if (!file_name_value.empty()) {
            // don't complain about defined_in rules
            continue;
         }

         const char* attrName = "class";
         const char* attrVal = nullptr;
         if (!name.empty()) attrVal = name.c_str();

         ROOT::TMetaUtils::Warning(0,"Unused %s rule: %s\n", attrName, attrVal);
      }
   }

   for(auto&& rule : fVariableSelectionRules){
      if (!rule.GetMatchFound() && !GetHasFileNameRule()) {
         std::string name;
         if (rule.GetAttributeValue("pattern", name)) {
            // keep it
         } else if (rule.GetAttributeValue("name", name)) {
            // keept it
         } else {
            name.clear();
         }
         ROOT::TMetaUtils::Warning("","Unused variable rule: %s\n",name.c_str());
         if (name.empty()) {
            rule.PrintAttributes(std::cout,3);
         }
      }
   }

#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
ROOT::TMetaUtils::Warning("SelectionRules::AreAllSelectionRulesUsed",
"Warnings concerning non matching selection rules are suppressed. An action is to be taken.\n");
#endif
#endif
//    for(const auto& selRule: fFunctionSelectionRules) {
//       if (!selRule.GetMatchFound() && !GetHasFileNameRule()) {
//          // Here the slow methods can be used
//          std::string name;
//          if (selRule.GetAttributeValue("proto_pattern", name)) {
//             // keep it
//          } else if (selRule.GetAttributeValue("proto_name", name)) {
//             // keep it
//          } else if (selRule.GetAttributeValue("pattern", name)) {
//             // keep it
//          } else if (selRule.GetAttributeValue("name", name)) {
//             // keept it
//          } else {
//             name.clear();
//          }
//          // Make it soft, no error - just warnings
//          std::cout<<"Warning - unused function rule: "<<name<<std::endl;
// //          if (IsSelectionXMLFile()){
// //             std::cout<<"Warning - unused function rule: "<<name<<std::endl;
// //          }
// //          else {
// //             std::cout<<"Error - unused function rule: "<<name<<std::endl;
// //          }
//          if (name.length() == 0) {
//             selRule.PrintAttributes(std::cout,3);
//          }
//       }
//
//    }


#if Enums_rules_becomes_useful_for_rootcling
   for(auto&& rule : fEnumSelectionRules){
      if (!rule.GetMatchFound() && !GetHasFileNameRule()) {
         std::string name;
         if (rule.GetAttributeValue("pattern", name)) {
            // keep it
         } else if (rule.GetAttributeValue("name", name)) {
            // keept it
         } else {
            name.clear();
         }

         ROOT::TMetaUtils::Warning("","Unused enum rule: %s\n",name.c_str());

         if (name.empty()){
            rule.PrintAttributes(std::cout,3);
         }
      }
   }
#endif
   return true;
}

bool SelectionRules::SearchNames(cling::Interpreter &interp)
{
   // std::cout<<"Searching Names In Selection Rules:"<<std::endl;
   for(std::list<ClassSelectionRule>::iterator it = fClassSelectionRules.begin(),
          end = fClassSelectionRules.end();
       it != end;
       ++it) {
      if (it->HasAttributeWithName("name")) {
         std::string name_value;
         it->GetAttributeValue("name", name_value);
         // In Class selection rules, we should be interested in scopes.
         const clang::Type *typeptr = 0;
         const clang::CXXRecordDecl *target
            = ROOT::TMetaUtils::ScopeSearch(name_value.c_str(), interp,
                                            true /*diag*/, &typeptr);
         if (target) {
            it->SetCXXRecordDecl(target,typeptr);
         }
      }
   }
   return true;
}


void SelectionRules::FillCache()
{
   // Fill the cache of every selection rule
   for (auto& rule : fClassSelectionRules) rule.FillCache();
   for (auto& rule : fFunctionSelectionRules) rule.FillCache();
   for (auto& rule : fVariableSelectionRules) rule.FillCache();
   for (auto& rule : fEnumSelectionRules) rule.FillCache();
}


