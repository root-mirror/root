// @(#)root/core/utils:$Id: SelectionRules.h 28529 2009-05-11 16:43:35Z pcanal $
// Author: Velislava Spasova September 2010

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef R__SELECTIONRULES_H
#define R__SELECTIONRULES_H

#include <list>
#include <string>
#include <vector>
#include <utility>

#include "BaseSelectionRule.h"
#include "ClassSelectionRule.h"
#include "VariableSelectionRule.h"
#include "clang/AST/Decl.h"

#include "TClingUtils.h"

namespace cling {
   class Interpreter;
}

namespace ROOT{
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}
#include <iostream>
namespace SelectionRulesUtils {

   template<class ASSOCIATIVECONTAINER>
   inline bool areEqualAttributes(const ASSOCIATIVECONTAINER& c1, const ASSOCIATIVECONTAINER& c2, bool moduloNameOrPattern){
      if (c1.size() != c2.size()) return false;
      if (moduloNameOrPattern) {
         for (auto&& keyValPairC1 : c1){
            auto keyC1 = keyValPairC1.first;
            if ("pattern" == keyC1 || "name" == keyC1) continue;
            auto valC1 = keyValPairC1.second;
            auto C2It = c2.find(keyC1);
            if (C2It == c2.end() || valC1 != C2It->second) return false;
         }
      }
      else {
         return !(c1 != c2);
      }
      return true;
   }

   template<class RULE>
   inline bool areEqual(const RULE* r1, const RULE* r2, bool moduloNameOrPattern = false){
      return areEqualAttributes(r1->GetAttributes(), r2->GetAttributes(), moduloNameOrPattern);
   }

   template<class RULESCOLLECTION>
   inline bool areEqualColl(const RULESCOLLECTION& r1,
                            const RULESCOLLECTION& r2,
                            bool moduloNameOrPattern = false){
      if (r1.size() != r2.size()) return false;
      auto rIt1 = r1.begin();
      auto rIt2 = r2.begin();
      for (;rIt1!=r1.cend();++rIt1,++rIt2){
         if (!areEqual(&(*rIt1),&(*rIt2), moduloNameOrPattern)) return false;
      }
      return true;
   }
   template<>
   inline bool areEqual<ClassSelectionRule>(const ClassSelectionRule* r1,
                                            const ClassSelectionRule* r2,
                                            bool moduloNameOrPattern){
      if (!areEqualAttributes(r1->GetAttributes(), r2->GetAttributes(),moduloNameOrPattern)) return false;
      // Now check fields
      if (!areEqualColl(r1->GetFieldSelectionRules(),
                        r2->GetFieldSelectionRules(),
                        true)) return false;
      // On the same footing, now check methods
      if (!areEqualColl(r1->GetMethodSelectionRules(),
                        r2->GetMethodSelectionRules(),
                        true)) return false;
      return true;
   }
}


class SelectionRules {

public:
   /// Type of selection file
   enum ESelectionFileTypes {
      kSelectionXMLFile,
      kLinkdefFile,
      kNumSelectionFileTypes
   };

   SelectionRules(cling::Interpreter &interp,
                  ROOT::TMetaUtils::TNormalizedCtxt& normCtxt,
                  const std::vector<std::pair<std::string,std::string>>& namesForExclusion):
      fSelectionFileType(kNumSelectionFileTypes),
      fHasFileNameRule(false),
      fRulesCounter(0),
      fNormCtxt(normCtxt),
      fInterp(interp) {
         long counter=1;
         for (auto& attrValPair : namesForExclusion){
            ClassSelectionRule csr(counter++, fInterp);
            csr.SetAttributeValue(attrValPair.first, attrValPair.second);
            csr.SetSelected(BaseSelectionRule::kNo);
            AddClassSelectionRule(csr);
            }
      }

   void AddClassSelectionRule(const ClassSelectionRule& classSel);
   bool HasClassSelectionRules() const { return !fClassSelectionRules.empty(); }
   const std::list<ClassSelectionRule>& GetClassSelectionRules() const {
      return fClassSelectionRules;
   }

   void AddFunctionSelectionRule(const FunctionSelectionRule& funcSel);
   bool HasFunctionSelectionRules() const {
      return !fFunctionSelectionRules.empty();
   }
   const std::list<FunctionSelectionRule>& GetFunctionSelectionRules() const {
      return fFunctionSelectionRules;
   }

   void AddVariableSelectionRule(const VariableSelectionRule& varSel);

   bool HasVariableSelectionRules() const {
      return !fVariableSelectionRules.empty();
   }
   const std::list<VariableSelectionRule>& GetVariableSelectionRules() const {
      return fVariableSelectionRules;
   }

   void AddEnumSelectionRule(const EnumSelectionRule& enumSel);
   bool HasEnumSelectionRules() const { return !fEnumSelectionRules.empty(); }
   const std::list<EnumSelectionRule>& GetEnumSelectionRules() const {
      return fEnumSelectionRules;
   }

   void PrintSelectionRules() const; // print all selection rules

   void ClearSelectionRules(); // clear all selection rules

   void SetHasFileNameRule(bool file_rule) { fHasFileNameRule = file_rule; }
   bool GetHasFileNameRule() const { return fHasFileNameRule; }

   int CheckDuplicates();
   void Optimize();

   // These method are called from clr-scan and return true if the Decl selected, false otherwise
   //const BaseSelectionRule  *IsDeclSelected(clang::Decl* D) const;
   const ClassSelectionRule *IsDeclSelected(const clang::RecordDecl* D, bool includeTypedefRule) const;
   const ClassSelectionRule *IsDeclSelected(const clang::TypedefNameDecl* D) const;
   const ClassSelectionRule *IsDeclSelected(const clang::NamespaceDecl* D) const;
   const BaseSelectionRule *IsDeclSelected(const clang::EnumDecl* D) const;
   const BaseSelectionRule *IsDeclSelected(const clang::VarDecl* D) const;
   const BaseSelectionRule *IsDeclSelected(const clang::FieldDecl* D) const;
   const BaseSelectionRule *IsDeclSelected(const clang::FunctionDecl* D) const;
   const BaseSelectionRule *IsDeclSelected(const clang::Decl* D) const;

   const ClassSelectionRule *IsClassSelected(const clang::Decl* D, const std::string& qual_name, bool includeTypedefRule) const; // is the class selected
   const ClassSelectionRule *IsNamespaceSelected(const clang::Decl* D, const std::string& qual_name) const; // is the class selected

   // is the global function, variable, enum selected - the behavior is different for linkdef.h and selection.xml - that's why
   // we have two functions
   const BaseSelectionRule *IsVarSelected(const clang::VarDecl* D, const std::string& qual_name) const;
   const BaseSelectionRule *IsFunSelected(const clang::FunctionDecl* D, const std::string& qual_name) const;
   const BaseSelectionRule *IsEnumSelected(const clang::EnumDecl* D, const std::string& qual_name) const;
   const BaseSelectionRule *IsLinkdefVarSelected(const clang::VarDecl* D, const std::string& qual_name) const;
   const BaseSelectionRule *IsLinkdefFunSelected(const clang::FunctionDecl* D, const std::string& qual_name) const;
   const BaseSelectionRule *IsLinkdefEnumSelected(const clang::EnumDecl* D, const std::string& qual_name) const;

   // is member (field, method, enum) selected; the behavior for linkdef.h methods is different
   const BaseSelectionRule *IsMemberSelected(const clang::Decl* D, const std::string& str_name) const;
   const BaseSelectionRule *IsLinkdefMethodSelected(const clang::Decl* D, const std::string& qual_name) const;

   // Return the number of rules
    unsigned int Size() const{return fClassSelectionRules.size()+
                                     fFunctionSelectionRules.size()+
                                     fVariableSelectionRules.size()+
                                     fEnumSelectionRules.size();};

   // returns true if the parent is class or struct
   bool IsParentClass(const clang::Decl* D) const;

   // the same but returns also the parent name and qualified name
   bool IsParentClass(const clang::Decl* D, std::string& parent_name, std::string& parent_qual_name) const;

   // returns the parent name and qualified name
   bool GetParentName(const clang::Decl* D, std::string& parent_name, std::string& parent_qual_name) const;


   //bool getParent(clang::Decl* D, clang::Decl* parent); - this method would have saved a lot of efforts but it crashes
   // and I didn't understand why

   // gets the name and qualified name of the Decl
   bool GetDeclName(const clang::Decl* D, std::string& name, std::string& qual_name) const;

   // gets the qualname of the decl, no checks performed
   void GetDeclQualName(const clang::Decl* D, std::string& qual_name) const;

   // gets the function prototype if the Decl (if it is global function or method)
   bool GetFunctionPrototype(const clang::FunctionDecl* F, std::string& prototype) const;

   bool IsSelectionXMLFile() const {
      return fSelectionFileType == kSelectionXMLFile;
   }
   bool IsLinkdefFile() const {
      return fSelectionFileType == kLinkdefFile;
   }
   void SetSelectionFileType(ESelectionFileTypes fileType) {
      fSelectionFileType = fileType;
   }

   // returns true if all selection rules are used at least once
   bool AreAllSelectionRulesUsed() const;

   // Go through all the selections rules and lookup the name if any in the AST.
   // and force the instantiation of template if any are used in the rules.
   bool SearchNames(cling::Interpreter &interp);

   void FillCache(); // Fill the cache of all selection rules

private:
   std::list<ClassSelectionRule>    fClassSelectionRules;    ///< List of the class selection rules
   std::list<FunctionSelectionRule> fFunctionSelectionRules; ///< List of the global functions selection rules
   std::list<VariableSelectionRule> fVariableSelectionRules; ///< List of the global variables selection rules
   std::list<EnumSelectionRule>     fEnumSelectionRules;     ///< List of the enums selection rules

   ESelectionFileTypes fSelectionFileType;

   bool fHasFileNameRule; ///< if we have a file name rule, this should be set to true
   long int fRulesCounter;

   ROOT::TMetaUtils::TNormalizedCtxt& fNormCtxt;
   cling::Interpreter &fInterp;

};

#endif
