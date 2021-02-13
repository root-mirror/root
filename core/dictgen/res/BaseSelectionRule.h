// @(#)root/core/utils:$Id: BaseSelectionRule.h 28529 2009-05-11 16:43:35Z pcanal $
// Author: Velislava Spasova September 2010

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef R__BASESELECTIONRULE_H
#define R__BASESELECTIONRULE_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BaseSelectionRule                                                    //
//                                                                      //
// Base selection class from which all                                  //
// selection classes should be derived                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <unordered_map>
#include <list>
#include <iosfwd>

namespace clang {
   class NamedDecl;
   class CXXRecordDecl;
   class Type;
}
namespace cling {
   class Interpreter;
}

class BaseSelectionRule
{
public:
   typedef std::unordered_map<std::string, std::string> AttributesMap_t; // The liste of selection rule's attributes (name, pattern, ...)

   enum ESelect { // a rule could be selected, vetoed or we don't care about it
      kYes,
      kNo,
      kDontCare
   };
   enum EMatchType {
      kName,
      kPattern,
      kFile,
      kNoMatch
   };

private:
   long                   fIndex;                  // Index indicating the ordering of the rules.
   long                   fLineNumber=-1;          // Line number of the selection file where the rule is located
   std::string            fSelFileName="";         // Name of the selection file
   AttributesMap_t        fAttributes;             // list of the attributes of the selection/exclusion rule
   ESelect                fIsSelected;             // selected/vetoed/don't care
   std::list<std::string> fSubPatterns;            // a list of subpatterns, generated form a pattern/proto_pattern attribute
   std::list<std::string> fFileSubPatterns;        // a list of subpatterns, generated form a file_pattern attribute
   bool                   fMatchFound;             // this is true if this selection rule has been used at least once
   const clang::CXXRecordDecl  *fCXXRecordDecl;    // Record decl of the entity searched for.
   const clang::Type           *fRequestedType;    // Same as the record decl but with some of the typedef preserved (Double32_t, Float16_t, etc..)
   cling::Interpreter *fInterp;

   // Cached for performance
   std::string fName = "";
   std::string fPattern = "";
   std::string fProtoName = "";
   std::string fProtoPattern = "";
   std::string fFileName = "";
   std::string fFilePattern = "";
   std::string fNArgsToKeep = "";
   bool fHasNameAttribute = false;
   bool fHasProtoNameAttribute = false;
   bool fHasPatternAttribute = false;
   bool fHasProtoPatternAttribute = false;
   bool fHasFileNameAttribute = false;
   bool fHasFilePatternAttribute = false;
   bool fHasFromTypedefAttribute = false;
   bool fIsFromTypedef = false;

public:

   BaseSelectionRule(ESelect sel) : fIndex(-1),fIsSelected(sel),fMatchFound(false),fCXXRecordDecl(NULL),fRequestedType(NULL),fInterp(NULL) {}

   BaseSelectionRule(long index, cling::Interpreter &interp, const char* selFileName = "", long lineno=-1) : fIndex(index),fLineNumber(lineno),fSelFileName(selFileName),fIsSelected(kNo),fMatchFound(false),fCXXRecordDecl(0),fRequestedType(0),fInterp(&interp) {}

   BaseSelectionRule(long index, ESelect sel, const std::string& attributeName, const std::string& attributeValue, cling::Interpreter &interp, const char* selFileName = "",long lineno=-1);

   virtual ~BaseSelectionRule() = default;

   virtual void DebugPrint() const;
   virtual void Print(std::ostream &out) const = 0;

   long    GetIndex() const { return fIndex; }
   void    SetIndex(long index) { fIndex=index; }

   long    GetLineNumber() const { return fLineNumber; }
   const char* GetSelFileName() const { return fSelFileName.c_str(); }

   bool    HasAttributeWithName(const std::string& attributeName) const; // returns true if there is an attribute with the specified name

   void    FillCache(); // Fill the cache for performant attribute retrival

   bool    GetAttributeValue(const std::string& attributeName, std::string& returnValue) const; // returns the value of the attribute with name attributeName

   inline const std::string& GetAttributeName() const {return fName;};
   inline bool HasAttributeName() const {return fHasNameAttribute;};

   inline const std::string& GetAttributeProtoName() const {return fProtoName;};
   inline bool HasAttributeProtoName() const {return fHasProtoNameAttribute;};

   inline const std::string& GetAttributePattern() const {return fPattern;};
   inline bool HasAttributePattern() const {return fHasPatternAttribute;};

   inline const std::string& GetAttributeProtoPattern() const {return fProtoPattern;};
   inline bool HasAttributeProtoPattern() const {return fHasProtoPatternAttribute;};

   inline const std::string& GetAttributeFileName() const {return fFileName;};
   inline bool HasAttributeFileName() const {return fHasFileNameAttribute;};

   inline const std::string& GetAttributeFilePattern() const {return fFilePattern;};
   inline bool HasAttributeFilePattern() const {return fHasFilePatternAttribute;};

   inline bool IsFromTypedef() const {return fIsFromTypedef;};
   inline bool HasAttributeFromTypedef() const {return fHasFromTypedefAttribute;};

   inline const std::string& GetAttributeNArgsToKeep() const {return fNArgsToKeep;};

   void    SetAttributeValue(const std::string& attributeName, const std::string& attributeValue); // sets an attribute with name attribute name and value attributeValue

   ESelect GetSelected() const;
   void    SetSelected(ESelect sel);

   bool  HasInterpreter() const {return fInterp!=NULL; };
   void  SetInterpreter(cling::Interpreter& interp) {fInterp=&interp; };

   const AttributesMap_t& GetAttributes() const; // returns the list of attributes
   void  PrintAttributes(int level) const;       // prints the list of attributes - level is the number of tabs from the beginning of the line
   void  PrintAttributes(std::ostream &out, int level) const;       // prints the list of attributes - level is the number of tabs from the beginning of the line

   EMatchType Match(const clang::NamedDecl *decl, const std::string& name, const std::string& prototype, bool isLinkdef) const; // for more detailed description look at the .cxx file

   void  SetMatchFound(bool match); // set fMatchFound
   bool  GetMatchFound() const;     // get fMatchFound

   const clang::Type *GetRequestedType() const;
   inline const clang::CXXRecordDecl *GetCXXRecordDecl() const {return fCXXRecordDecl;} ;
   void SetCXXRecordDecl(const clang::CXXRecordDecl *decl, const clang::Type *typeptr);

protected:

   // Checks if the test string matches against the pattern (which has produced the list of sub-patterns patterns_list). There is
   // difference if we are processing linkdef.h or selection.xmlpatterns
   inline bool CheckPattern(const std::string& test, const std::string& pattern, const std::list<std::string>& patterns_list, bool isLinkdef) const;

   inline void ProcessPattern(const std::string& pattern, std::list<std::string>& out) const; // divides a pattern into a list of sub-patterns
};

inline std::ostream &operator<<(std::ostream& out, const BaseSelectionRule &obj)
{
   obj.Print(out);
   return out;
}
#endif
