// @(#)root/core/utils:$Id: BaseSelectionRule.cxx 41697 2011-11-01 21:03:41Z pcanal $
// Author: Velislava Spasova September 2010

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BaseSelectionRule                                                    //
//                                                                      //
// Base selection class from which all                                  //
// selection classes should be derived                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "BaseSelectionRule.h"

#include "TClassEdit.h"
#include "TClingUtils.h"

#include <iostream>
#include <string.h>
#include <cctype>

#include "clang/Basic/SourceLocation.h"
#include "clang/Basic/SourceManager.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/DeclTemplate.h"

#ifdef _WIN32
#include "process.h"
#endif
#include <sys/stat.h>

static const char *R__GetDeclSourceFileName(const clang::Decl* D)
{
   clang::ASTContext& ctx = D->getASTContext();
   clang::SourceManager& SM = ctx.getSourceManager();
   clang::SourceLocation SL = D->getLocation();
   // If the class decl is the result of a macpo expansion, take the location
   // where the macro is "invoked" i.e. expanded at (ExpansionLoc), not the
   // spelling location (where the delc's tokens come from).
   if (SL.isMacroID())
      SL = SM.getExpansionLoc(SL);

   if (SL.isValid() && SL.isFileID()) {
      clang::PresumedLoc PLoc = SM.getPresumedLoc(SL);
      return PLoc.getFilename();
   }
   else {
      return "invalid";
   }
}

static bool R__match_filename(const char *srcname,const char *filename)
{
   if (srcname==0) {
      return false;
   }
   if((strcmp(srcname,filename)==0)) {
      return true;
   }

#ifdef G__WIN32
   char i1name[_MAX_PATH];
   char fullfile[_MAX_PATH];
   _fullpath( i1name, srcname, _MAX_PATH );
   _fullpath( fullfile, filename, _MAX_PATH );
   if((stricmp(i1name, fullfile)==0)) return 1;
#else
   struct stat statBufItem;
   struct stat statBuf;
   if (   ( 0 == stat( filename, & statBufItem ) )
       && ( 0 == stat( srcname, & statBuf ) )
       && ( statBufItem.st_dev == statBuf.st_dev )     // Files on same device
       && ( statBufItem.st_ino == statBuf.st_ino )     // Files on same inode (but this is not unique on AFS so we need the next 2 test
       && ( statBufItem.st_size == statBuf.st_size )   // Files of same size
       && ( statBufItem.st_mtime == statBuf.st_mtime ) // Files modified at the same time
       ) {
      return true;
   }
#endif
   return false;
}

BaseSelectionRule::BaseSelectionRule(long index, BaseSelectionRule::ESelect sel, const std::string& attributeName, const std::string& attributeValue, cling::Interpreter &interp, const char* selFileName, long lineno)
   : fIndex(index),fLineNumber(lineno),fSelFileName(selFileName),fIsSelected(sel),fMatchFound(false),fCXXRecordDecl(0),fRequestedType(0),fInterp(&interp)
{
   fAttributes.insert(AttributesMap_t::value_type(attributeName, attributeValue));
}

void BaseSelectionRule::SetSelected(BaseSelectionRule::ESelect sel)
{
   fIsSelected = sel;
}

BaseSelectionRule::ESelect BaseSelectionRule::GetSelected() const
{
   return fIsSelected;
}

bool BaseSelectionRule::HasAttributeWithName(const std::string& attributeName) const
{
   AttributesMap_t::const_iterator iter = fAttributes.find(attributeName);

   if(iter!=fAttributes.end()) return true;
   else return false;
}

bool BaseSelectionRule::GetAttributeValue(const std::string& attributeName, std::string& returnValue) const
{
   AttributesMap_t::const_iterator iter = fAttributes.find(attributeName);

   bool retVal = iter!=fAttributes.end();
   returnValue = retVal ? iter->second : "";
   return retVal;
}

void BaseSelectionRule::SetAttributeValue(const std::string& attributeName, const std::string& attributeValue)
{

   std::string localAttributeValue(attributeValue);

   int pos = attributeName.find("pattern");
   int pos_file = attributeName.find("file_pattern");

   // Strip trailing spaces from the name or pattern
   if (attributeName == "name" || pos> -1){
      while(std::isspace(*localAttributeValue.begin())) localAttributeValue.erase(localAttributeValue.begin());
      while(std::isspace(*localAttributeValue.rbegin()))localAttributeValue.erase(localAttributeValue.length()-1);
   }
   fAttributes.insert(AttributesMap_t::value_type(attributeName, localAttributeValue));

   if (pos > -1) {
      if (pos_file > -1) // if we have file_pattern
         ProcessPattern(localAttributeValue, fFileSubPatterns);
      else ProcessPattern(localAttributeValue, fSubPatterns); // if we have pattern and proto_pattern
   }



}

const BaseSelectionRule::AttributesMap_t& BaseSelectionRule::GetAttributes() const
{
   return fAttributes;
}

void BaseSelectionRule::DebugPrint() const
{
   Print(std::cout);
}

void BaseSelectionRule::PrintAttributes(std::ostream &out, int level) const
{
   std::string tabs;
   for (int i = 0; i < level; ++i) {
      tabs+='\t';
   }

   if (!fAttributes.empty()) {
      std::map<std::string,std::string> orderedAttributes(fAttributes.begin(),fAttributes.end());
      for (auto&& attr : orderedAttributes) {
         out<<tabs<<attr.first<<" = "<<attr.second<<std::endl;
      }
   }
   else {
      out<<tabs<<"No attributes"<<std::endl;
   }
}

void BaseSelectionRule::PrintAttributes(int level) const
{
   PrintAttributes(std::cout, level);
}
#ifndef G__WIN32
#include <unistd.h>
#endif
BaseSelectionRule::EMatchType BaseSelectionRule::Match(const clang::NamedDecl *decl,
                                                       const std::string& name,
                                                       const std::string& prototype,
                                                       bool isLinkdef) const
{
   /* This method returns whether and how the declaration is matching the rule.
    * It returns one of:
    *   kNoMatch : the rule does match the declaration
    *   kName    : the rule match the declaration by name
    *   kPattern : the rule match the declaration via a pattern
    *   kFile    : the declaration's file name is match by the rule (either by name or pattern).
    * To check whether the rule is accepting or vetoing the declaration see the result of
    * GetSelected().
    (
    * We pass as arguments of the method:
    *   name - the name of the Decl
    *   prototype - the prototype of the Decl (if it is function or method, otherwise "")
    *   file_name - name of the source file
    *   isLinkdef - if the selection rules were generating from a linkdef.h file
    */

   const std::string& name_value = fName;
   const std::string& pattern_value = fPattern;

   // Check if we have in hands a typedef to a RecordDecl
   const clang::CXXRecordDecl *D = llvm::dyn_cast<clang::CXXRecordDecl>(decl);
   bool isTypedefNametoRecordDecl = false;

   if (!D){
      //Either it's a CXXRecordDecl ot a TypedefNameDecl
      const clang::TypedefNameDecl* typedefNameDecl = llvm::dyn_cast<clang::TypedefNameDecl> (decl);
      isTypedefNametoRecordDecl = typedefNameDecl &&
                                  ROOT::TMetaUtils::GetUnderlyingRecordDecl(typedefNameDecl->getUnderlyingType());
      }

   if (! isTypedefNametoRecordDecl && fCXXRecordDecl !=0 && fCXXRecordDecl != (void*)-1) {
      const clang::CXXRecordDecl *target = fCXXRecordDecl;
      if ( target && D && target == D ) {
         //               fprintf(stderr,"DECL MATCH: %s %s\n",name_value.c_str(),name.c_str());
         const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
         return kName;
      }
   } else if (fHasNameAttribute) {
      if (name_value == name) {
         const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
         return kName;
      } else if (fCXXRecordDecl == nullptr ||
            (fCXXRecordDecl != (void*)-1 && isTypedefNametoRecordDecl && !decl->hasOwningModule())){
         // Possibly take the most expensive path if the fCXXRecordDecl is not
         // set or we already took the expensive path and found nothing (-1).
         const clang::CXXRecordDecl *target
            = fHasFromTypedefAttribute ? nullptr : ROOT::TMetaUtils::ScopeSearch(name_value.c_str(), *fInterp,
                                                   true /*diagnose*/, 0);

         if ( target ) {
            const_cast<BaseSelectionRule*>(this)->fCXXRecordDecl = target;
         } else {
            // If the lookup failed, let's not try it again, so mark the value has invalid.
            const_cast<BaseSelectionRule*>(this)->fCXXRecordDecl = (clang::CXXRecordDecl*)-1;
         }
         if ( target && D && target == D ) {
            const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
            return kName;
         }
      }
   }

   // do we have matching against the file_name (or file_pattern) attribute and if yes - select or veto
   const std::string& file_name_value = fFileName;
   const std::string& file_pattern_value = fFilePattern;

   if ((fHasFileNameAttribute||fHasFilePatternAttribute)) {
      const char *file_name = R__GetDeclSourceFileName(decl);
      bool hasFileMatch = ((fHasFileNameAttribute &&
           //FIXME It would be much better to cache the rule stat result and compare to the clang::FileEntry
           (R__match_filename(file_name_value.c_str(),file_name))) ||
           (fHasFilePatternAttribute && CheckPattern(file_name, file_pattern_value, fFileSubPatterns, isLinkdef)));

      if (hasFileMatch) {
         // Reject utility classes defined in ClassImp
         // when using a file based rule
         if (!strncmp(name.c_str(), "R__Init", 7) ||
             strstr(name.c_str(), "::R__Init")) {
            return kNoMatch;
         }
         if (!name.compare(0, 24, "ROOT::R__dummyintdefault")) {
            return kNoMatch;
         }
         if (!name.compare(0, 27, "ROOT::R__dummyVersionNumber")) {
            return kNoMatch;
         }
         if (!name.compare(0, 22, "ROOT::R__dummyStreamer")) {
            return kNoMatch;
         }
         if (name.find("(anonymous namespace)") != std::string::npos) {
            // Reject items declared in anonymous namespace
            return kNoMatch;
         }
         if (fHasPatternAttribute) {
            if (CheckPattern(name, pattern_value, fSubPatterns, isLinkdef)) {
               const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
               return kPattern;
            }
         } else {
            const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
            return kName;
         }
      }

      // We have file_name or file_pattern attribute but the
      // passed file_name is different than that in the selection rule then return no match
      return kNoMatch;
   }

   if (fHasPatternAttribute)
   {
      bool patternMatched = CheckPattern(name, pattern_value, fSubPatterns, isLinkdef);
      if (!patternMatched && !isLinkdef) {
         std::string auxName(name);
         std::string &nameNoSpaces = auxName;
         nameNoSpaces.erase(std::remove_if(nameNoSpaces.begin(), nameNoSpaces.end(), isspace),
                           nameNoSpaces.end());
         if (name.size() != nameNoSpaces.size()) {
            patternMatched = CheckPattern(nameNoSpaces, pattern_value, fSubPatterns, isLinkdef);
         }

         // For ROOT-6704: use normalised name for matching if the class is in stl
         // The reason for this check is that we have rules like std::map<*, int>
         // We do not know how the internal representation of the innocuous "map"
         // is. We therefore have to act on a nicer name, obtained with TClassEdit
         // The check ROOT::TMetaUtils::IsStdDropDefaultClass is there to call
         // TClassEdit only when necessary as it can be expensive, a performance
         // optimisation.
         if (!patternMatched &&
               D &&
               //ROOT::TMetaUtils::IsStdDropDefaultClass(*D)) {
               0 != TClassEdit::IsSTLCont(name)) {
            TClassEdit::GetNormalizedName(auxName, name);
            if (name.size() != auxName.size()) {
               auxName = TClassEdit::InsertStd(auxName.c_str());
               patternMatched = CheckPattern(auxName, pattern_value, fSubPatterns, isLinkdef);
            }
         }

      }
      if (patternMatched) {
         const_cast<BaseSelectionRule *>(this)->SetMatchFound(true);
         return kPattern;
      }
   }


   // do we have matching against the proto_name (or proto_pattern)  attribute and if yes - select or veto
   // The following selects functions on whether the requested prototype exactly matches the
   // prototype issued by SelectionRules::GetFunctionPrototype which relies on
   //    ParmVarDecl::getType()->getAsString()
   // to get the type names.  Currently, this does not print the prototype in the usual
   // human (written) forms.   For example:
   //   For Hash have prototype: '(const class TString &)'
   //   For Hash have prototype: '(const class TString*)'
   //   For Hash have prototype: '(const char*)'
   // In addition, the const can legally be in various place in the type name and thus
   // a string based match will be hard to work out (it would need to normalize both
   // the user input string and the clang provided string).
   // Using lookup form cling would be probably be a better choice.
   if (!prototype.empty()) {
      if (fHasProtoNameAttribute && fProtoName==prototype) {
         const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
         return kName;
      }
      if (fHasProtoPatternAttribute && CheckPattern(prototype, fProtoPattern, fSubPatterns, isLinkdef))  {
         const_cast<BaseSelectionRule*>(this)->SetMatchFound(true);
         return kPattern;
      }
   }

   return kNoMatch;
}


/*
 * This method processes the pattern - which means that it splits it in a list of fSubPatterns.
 * The idea is the following - if we have a pattern = "this*pat*rn", it will be split in the
 * following list of subpatterns: "this", "pat", "rn". If we have "this*pat\*rn", it will be
 * split in "this", "pat*rn", i.e. the star could be escaped.
 */

void BaseSelectionRule::ProcessPattern(const std::string& pattern, std::list<std::string>& out) const
{
   std::string temp = pattern;
   std::string split;
   int pos;
   bool escape = false;

   if (pattern.size()==1 && pattern == "*"){
      out.push_back("");
      return;
   }

   while (!temp.empty()){
      pos = temp.find("*");
      if (pos == -1) {
         if (!escape){ // if we don't find a '*', push_back temp (contains the last sub-pattern)
            out.push_back(temp);
            // std::cout<<"1. pushed = "<<temp<<std::endl;
         }
         else { // if we don't find a star - add temp to split (in split we keep the previous sub-pattern + the last escaped '*')
            split += temp;
            out.push_back(split);
            // std::cout<<"1. pushed = "<<split<<std::endl;
         }
         return;
      }
      else if (pos == 0) { // we have '*' at the beginning of the pattern; can't have '\' before the '*'
         temp = temp.substr(1); // remove the '*'
      }
      else if (pos == (int)(temp.length()-1)) { // we have '*' at the end of the pattern
         if (pos > 0 && temp.at(pos-1) == '\\') { // check if we have '\' before the '*'; if yes, we have to escape it
            split += temp.substr(0, temp.length()-2);  // add evrything from the beginning of temp till the '\' to split (where we keep the last sub-pattern)
            split += temp.at(pos); // add the '*'
            out.push_back(split);  // push_back() split
            // std::cout<<"3. pushed = "<<split<<std::endl;
            temp.clear(); // empty temp (the '*' was at the last position of temp, so we don't have anything else to process)
         }
         temp = temp.substr(0, (temp.length()-1));
      }
      else { // the '*' is at a random position in the pattern
         if (pos > 0 && temp.at(pos-1) == '\\') { // check if we have '\' before the '*'; if yes, we have to escape it
            split += temp.substr(0, pos-1); // remove the '\' and add the star to split
            split += temp.at(pos);
            escape = true;                  // escape = true which means that we will add the next sub-pattern to that one

            // DEBUG std::cout<<"temp = "<<temp<<std::endl;
            temp = temp.substr(pos);
            // DEBUG std::cout<<"temp = "<<temp<<", split = "<<split<<std::endl;
         }
         else { // if we don't have '\' before the '*'
            if (escape) {
               split += temp.substr(0, pos);
            }
            else {
               split = temp.substr(0, pos);
            }
            escape = false;
            temp = temp.substr(pos);
            out.push_back(split);
            // std::cout<<"2. pushed = "<<split<<std::endl;
            // DEBUG std::cout<<"temp = "<<temp<<std::endl;
            split = "";
         }
      }
      // DEBUG std::cout<<"temp = "<<temp<<std::endl;
   }
}

/*
 * This method checks if the given test string is matched against the pattern
 */

bool BaseSelectionRule::CheckPattern(const std::string& test, const std::string& pattern, const std::list<std::string>& patterns_list, bool isLinkdef) const
{
   bool begin = pattern.front() == '*';
   if (pattern.size() == 1 && begin) {
      // We have the simple pattern '*', it matches everything by definition!
      return true;
   }

   std::list<std::string>::const_iterator it = patterns_list.begin();
   size_t pos1, pos2, pos3;
   pos1= pos2= pos3= std::string::npos;
   bool end = pattern.back() == '*';

   // we first check if the last sub-pattern is contained in the test string
   const std::string& last = patterns_list.back();
   size_t pos_end = test.rfind(last);

   if (pos_end == std::string::npos) { // the last sub-pattern isn't conatained in the test string
      return false;
   }
   if (!end) {  // if the pattern doesn't end with '*', the match has to be complete
      // i.e. if the last sub-pattern is "sub" the test string should end in "sub" ("1111sub" is OK, "1111sub1" is not OK)

      int len = last.length(); // length of last sub-pattern
      if ((pos_end+len) < test.length()) {
         return false;
      }
   }

   // position of the first sub-pattern
   pos1 = test.find(*it);


   if (pos1 == std::string::npos || (!begin && pos1 != 0)) { // if the first sub-pattern isn't found in test or if it is found but the
      // pattern doesn't start with '*' and the sub-pattern is not at the first position
      //std::cout<<"\tNo match!"<<std::endl;
      return false;
   }

   if (isLinkdef) { // A* selects all global classes, unions, structs but not the nested, i.e. not A::B
      // A::* selects the nested classes
      int len = (*it).length();
      int pos_colon = test.find("::", pos1+len);

      if (pos_colon > -1) {
         return false;
      }

   }

   if (patterns_list.size() > 1) {
      if (((*it).length())+pos1 > pos_end) {
         // std::cout<<"\tNo match";
         return false; // end is contained in begin -> test = "A::B" sub-patterns = "A::", "::" will return false
      }
   }


   ++it;

   for (; it != patterns_list.end(); ++it) {
      // std::cout<<"sub-pattern = "<<*it<<std::endl;
      pos2 = test.find(*it);
      if (pos2 <= pos1) {
         return false;
      }
      pos1 = pos2;
   }

   return true;
}


void BaseSelectionRule::SetMatchFound(bool match)
{
   fMatchFound = match;
}

bool BaseSelectionRule::GetMatchFound() const
{
   return fMatchFound;
}

const clang::Type *BaseSelectionRule::GetRequestedType() const
{
   return fRequestedType;
}

void BaseSelectionRule::SetCXXRecordDecl(const clang::CXXRecordDecl *decl, const clang::Type *typeptr)
{
   fCXXRecordDecl = decl;
   fRequestedType = typeptr;
}

void BaseSelectionRule::FillCache()
{
   std::string value;
   fHasNameAttribute = GetAttributeValue("name",fName);
   fHasProtoNameAttribute = GetAttributeValue("proto_name",fProtoName);
   fHasPatternAttribute = GetAttributeValue("pattern",fPattern);
   fHasProtoPatternAttribute = GetAttributeValue("proto_pattern",fProtoPattern);
   fHasFileNameAttribute = GetAttributeValue("file_name",fFileName);
   fHasFilePatternAttribute = GetAttributeValue("file_pattern",fFilePattern);
   fHasFromTypedefAttribute = GetAttributeValue("fromTypedef",value);
   fIsFromTypedef = (value == "true");

   GetAttributeValue(ROOT::TMetaUtils::propNames::nArgsToKeep,fNArgsToKeep);


   if (fHasPatternAttribute || fHasProtoPatternAttribute) {
      if (fSubPatterns.empty()) {
         std::cout<<"Error - A pattern selection without sub patterns." <<std::endl;
      }
   }

}


