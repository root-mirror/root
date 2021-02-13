// @(#)root/html:$Id$
// Author: Axel Naumann 2007-01-09

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TDocParser.h"

#include "Riostream.h"
#include "TBaseClass.h"
#include "TClass.h"
#include "TClassDocOutput.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TDatime.h"
#include "TDocDirective.h"
#include "TGlobal.h"
#include "THtml.h"
#include "TInterpreter.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TVirtualMutex.h"
#include <string>

namespace {

   class TMethodWrapperImpl: public TDocMethodWrapper {
   public:
      TMethodWrapperImpl(TMethod* m, int overloadIdx):
         fMeth(m), fOverloadIdx(overloadIdx) {}

      static void SetClass(const TClass* cl) { fgClass = cl; }

      const char* GetName() const { return fMeth->GetName(); }
      ULong_t Hash() const { return fMeth->Hash();}
      Int_t GetNargs() const { return fMeth->GetNargs(); }
      virtual TMethod* GetMethod() const { return fMeth; }
      Bool_t IsSortable() const { return kTRUE; }

      Int_t GetOverloadIdx() const { return fOverloadIdx; }

      Int_t Compare(const TObject *obj) const {
         const TMethodWrapperImpl* m = dynamic_cast<const TMethodWrapperImpl*>(obj);
         if (!m) return 1;

         Int_t ret = strcasecmp(GetName(), m->GetName());
         if (ret == 0) {
            if (GetNargs() < m->GetNargs()) return -1;
            else if (GetNargs() > m->GetNargs()) return 1;
            if (GetMethod()->GetClass()->InheritsFrom(m->GetMethod()->GetClass()))
               return -1;
            else
               return 1;
         }

         const char* l(GetName());
         const char* r(m->GetName());
         if (l[0] == '~' && r[0] == '~') {
            ++l;
            ++r;
         }
         TClass *lcl = 0;
         TClass *rcl = 0;
         if (fMeth->Property() & (kIsConstructor|kIsDestructor)) {
            lcl = TClass::GetClass(l);
         }
         if (m->fMeth->Property() & (kIsConstructor|kIsDestructor)) {
            rcl = TClass::GetClass(r);
         }
         if (lcl && fgClass->InheritsFrom(lcl)) {
               if (rcl && fgClass->InheritsFrom(rcl)) {
                  if (lcl->InheritsFrom(rcl))
                     return -1;
                  else return 1;
               } else return -1;
         } else if (rcl && fgClass->InheritsFrom(rcl))
            return 1;

         if (l[0] == '~') return -1;
         if (r[0] == '~') return 1;
         return (ret < 0) ? -1 : 1;
      }

   private:
      static const TClass* fgClass; // current class, defining inheritance sort order
      TMethod* fMeth; // my method
      Int_t fOverloadIdx; // this is the n-th overload
   };

   const TClass* TMethodWrapperImpl::fgClass = 0;
}


//______________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////
//
// Parse C++ source or header, and extract documentation.
//
// Also handles special macros like
/* Begin_Macro(GUI, source)
{
   TGMainFrame* f = new TGMainFrame(0, 100, 100);
   f->SetName("testMainFrame"); // that's part of the name of the image
   TGButton* b = new TGTextButton(f, "Test Button");
   f->AddFrame(b);
   f->MapSubwindows();
   f->Resize(f->GetDefaultSize());

   f->MapWindow();
   return f; // *HIDE*
}
End_Macro */
// or multiline Latex aligned at =:
/* Begin_Latex(separator='=',align=rcl) C = d #sqrt{#frac{2}{#lambdaD}} #int^{x}_{0}cos(#frac{#pi}{2}t^{2})dt
 D(x) = d End_Latex */
// even without alignment: Begin_Latex
// x=sin^2(y)
// y = #sqrt{sin(x)}
// End_Latex and what about running an external macro?
/* BEGIN_MACRO(source)


testmacro.C END_MACRO


and some nested stuff which doesn't work yet: */
// BEGIN_HTML
/* BEGIN_LATEX Wow,^{an}_{image}^{inside}_{a}^{html}_{block}
   END_LATEX
*/
// END_HTML
////////////////////////////////////////////////////////////////////////////////

ClassImp(TDocParser);

std::set<std::string>  TDocParser::fgKeywords;

////////////////////////////////////////////////////////////////////////////////
/// Constructor called for parsing class sources

TDocParser::TDocParser(TClassDocOutput& docOutput, TClass* cl):
   fHtml(docOutput.GetHtml()), fDocOutput(&docOutput), fLineNo(0),
   fCurrentClass(cl), fRecentClass(0), fCurrentModule(0),
   fDirectiveCount(0), fLineNumber(0), fDocContext(kIgnore),
   fCheckForMethod(kFALSE), fClassDocState(kClassDoc_Uninitialized),
   fCommentAtBOL(kFALSE), fAllowDirectives(kTRUE)
{
   InitKeywords();

   fSourceInfoTags[kInfoLastUpdate] = fHtml->GetLastUpdateTag();
   fSourceInfoTags[kInfoAuthor]     = fHtml->GetAuthorTag();
   fSourceInfoTags[kInfoCopyright]  = fHtml->GetCopyrightTag();

   fClassDescrTag = fHtml->GetClassDocTag();

   TMethodWrapperImpl::SetClass(cl);

   for (int ia = 0; ia < 3; ++ia) {
      fMethods[ia].Rehash(101);
   }

   AddClassMethodsRecursively(0);
   AddClassDataMembersRecursively(0);

   // needed for list of methods,...
   fParseContext.push_back(kCode);

   // create an array of method names
   TMethod *method;
   TIter nextMethod(fCurrentClass->GetListOfMethods());
   fMethodCounts.clear();
   while ((method = (TMethod *) nextMethod())) {
      ++fMethodCounts[method->GetName()];
   }

}

////////////////////////////////////////////////////////////////////////////////
/// constructor called for parsing text files with Convert()

TDocParser::TDocParser(TDocOutput& docOutput):
   fHtml(docOutput.GetHtml()), fDocOutput(&docOutput), fLineNo(0),
   fCurrentClass(0), fRecentClass(0), fDirectiveCount(0),
   fLineNumber(0), fDocContext(kIgnore),
   fCheckForMethod(kFALSE), fClassDocState(kClassDoc_Uninitialized),
   fCommentAtBOL(kFALSE), fAllowDirectives(kFALSE)
{
   InitKeywords();

   fSourceInfoTags[kInfoLastUpdate] = fHtml->GetLastUpdateTag();
   fSourceInfoTags[kInfoAuthor]     = fHtml->GetAuthorTag();
   fSourceInfoTags[kInfoCopyright]  = fHtml->GetCopyrightTag();

   fClassDescrTag = fHtml->GetClassDocTag();

   TMethodWrapperImpl::SetClass(0);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor, checking whether all methods have been found for gDebug > 3

TDocParser::~TDocParser()
{
   if (gDebug > 3) {
      for (std::map<std::string, Int_t>::const_iterator iMethod = fMethodCounts.begin();
         iMethod != fMethodCounts.end(); ++iMethod)
         if (iMethod->second)
            Info("~TDocParser", "Implementation of method %s::%s could not be found.",
            fCurrentClass ? fCurrentClass->GetName() : "",
            iMethod->first.c_str());
      TIter iDirective(&fDirectiveHandlers);
      TDocDirective* directive = 0;
      while ((directive = (TDocDirective*) iDirective())) {
         TString directiveName;
         directive->GetName(directiveName);
         Warning("~TDocParser", "Missing \"%s\" for macro %s", directive->GetEndTag(), directiveName.Data());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add accessible (i.e. non-private) methods of base class bc
/// and its base classes' methods to methodNames.
/// If bc==0, we add fCurrentClass's methods (and also private functions).

void TDocParser::AddClassMethodsRecursively(TBaseClass* bc)
{
   // make a loop on member functions
   TClass *cl = fCurrentClass;
   if (bc)
      cl = bc->GetClassPointer(kFALSE);
   if (!cl) return;

   TMethod *method;
   TIter nextMethod(cl->GetListOfMethods());
   std::map<std::string, int> methOverloads;

   while ((method = (TMethod *) nextMethod())) {

      if (!strcmp(method->GetName(), "Dictionary") ||
          !strcmp(method->GetName(), "Class_Version") ||
          !strcmp(method->GetName(), "Class_Name") ||
          !strcmp(method->GetName(), "DeclFileName") ||
          !strcmp(method->GetName(), "DeclFileLine") ||
          !strcmp(method->GetName(), "ImplFileName") ||
          !strcmp(method->GetName(), "ImplFileLine") ||
          (bc && (method->GetName()[0] == '~' // d'tor
             || !strcmp(method->GetName(), method->GetReturnTypeName()))) // c'tor
          )
         continue;


      Int_t mtype = 0;
      if (kIsPrivate & method->Property())
         mtype = 0;
      else if (kIsProtected & method->Property())
         mtype = 1;
      else if (kIsPublic & method->Property())
         mtype = 2;

      if (bc) {
         if (mtype == 0) continue;
         if (bc->Property() & kIsPrivate)
            mtype = 0;
         else if ((bc->Property() & kIsProtected) && mtype == 2)
            mtype = 1;
      }

      Bool_t hidden = kFALSE;
      for (Int_t access = 0; !hidden && access < 3; ++access) {
         TMethodWrapperImpl* other = (TMethodWrapperImpl*) fMethods[access].FindObject(method->GetName());
         hidden |= (other) && (other->GetMethod()->GetClass() != method->GetClass());
      }
      if (!hidden) {
         fMethods[mtype].Add(new TMethodWrapperImpl(method, methOverloads[method->GetName()]));
         ++methOverloads[method->GetName()];
      }
   }

   TIter iBase(cl->GetListOfBases());
   TBaseClass* base = 0;
   while ((base = (TBaseClass*)iBase()))
      AddClassMethodsRecursively(base);

   if (!bc)
      for (Int_t access = 0; access < 3; ++access) {
         fMethods[access].SetOwner();
         fMethods[access].Sort();
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Add data members of fCurrentClass and of bc to datamembers, recursively.
/// Real data members are in idx 0..2 (public, protected, private access),
/// enum constants in idx 3..5.

void TDocParser::AddClassDataMembersRecursively(TBaseClass* bc) {
   // make a loop on member functions
   TClass *cl = fCurrentClass;
   if (bc)
      cl = bc->GetClassPointer(kFALSE);
   if (!cl) return;

   TDataMember *dm;
   TIter nextDM(cl->GetListOfDataMembers());

   while ((dm = (TDataMember *) nextDM())) {
      if (!strcmp(dm->GetName(), "fgIsA"))
         continue;
      Int_t mtype = 0;
      if (kIsPrivate & dm->Property())
         mtype = 0;
      else if (kIsProtected & dm->Property())
         mtype = 1;
      else if (kIsPublic & dm->Property())
         mtype = 2;

      if (bc) {
         if (mtype == 0) continue;
         if (bc->Property() & kIsPrivate)
            mtype = 0;
         else if ((bc->Property() & kIsProtected) && mtype == 2)
            mtype = 1;
      }

      const Int_t flagEnumConst = kIsEnum | kIsConstant | kIsStatic;
      if ((dm->Property() & flagEnumConst) == flagEnumConst
          && dm->GetDataType() && dm->GetDataType()->GetType() == kInt_t) {
         mtype = 5;
         // The access of the enum constant is defined by the access of the enum:
         // for CINT, all enum constants are public.
         // There is no TClass or TDataType for enum types; instead, use CINT:
         /*
           No - CINT does not know their access restriction.
           With CINT5 we have no way of determining it...

         ClassInfo_t* enumCI = gInterpreter->ClassInfo_Factory(dm->GetTypeName());
         if (enumCI) {
            Long_t prop = gInterpreter->ClassInfo_Property(enumCI);
            if (kIsPrivate & prop)
               mtype = 3;
            else if (kIsProtected & prop)
               mtype = 4;
            else if (kIsPublic & prop)
               mtype = 5;
            gInterpreter->ClassInfo_Delete(enumCI);
         }
         */
      }

      fDataMembers[mtype].Add(dm);
   }

   TIter iBase(cl->GetListOfBases());
   TBaseClass* base = 0;
   while ((base = (TBaseClass*)iBase()))
      AddClassDataMembersRecursively(base);

   if (!bc)
      for (Int_t access = 0; access < 6; ++access) {
         fDataMembers[access].SetOwner(kFALSE);
         if (access < 3) // don't sort enums; we keep them in enum tag order
            fDataMembers[access].Sort();
      }
}


////////////////////////////////////////////////////////////////////////////////
/// Create an anchor from the given line, by hashing it and
/// convertig the hash into a custom base64 string.

void TDocParser::AnchorFromLine(const TString& line, TString& anchor) {
   const char base64String[65] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.";

   // use hash of line instead of e.g. line number.
   // advantages: more stable (lines can move around, we still find them back),
   // no need for keeping a line number context
   UInt_t hash = ::Hash(line);
   anchor.Remove(0);
   // force first letter to be [A-Za-z], to be id compatible
   anchor += base64String[hash % 52];
   hash /= 52;
   while (hash) {
      anchor += base64String[hash % 64];
      hash /= 64;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Parse text file "in", add links etc, and write output to "out".
/// If "isCode", "in" is assumed to be C++ code.

void TDocParser::Convert(std::ostream& out, std::istream& in, const char* relpath,
                         Bool_t isCode, Bool_t interpretDirectives)
{
   fLineNumber = 0;
   fParseContext.clear();
   if (isCode) fParseContext.push_back(kCode);
   else        fParseContext.push_back(kComment); // so we can find "BEGIN_HTML"/"END_HTML" in plain text

   while (!in.eof()) {
      fLineRaw.ReadLine(in, kFALSE);
      ++fLineNumber;
      if (in.eof())
         break;

      // remove leading spaces
      fLineComment = "";
      fLineSource = fLineRaw;
      fLineStripped = fLineRaw;
      Strip(fLineStripped);

      DecorateKeywords(fLineSource);

      // Changes in this bit of code have consequences for:
      // * module index,
      // * source files,
      // * THtml::Convert() e.g. in tutorials/html/MakeTutorials.C
      if (!interpretDirectives) {
         // Only write the raw, uninterpreted directive code:
         if (!InContext(kDirective)) {
            GetDocOutput()->AdjustSourcePath(fLineSource, relpath);
            out << fLineSource << std::endl;
         }
      } else {
         // Write source for source and interpreted directives if they exist.
         if (fLineComment.Length() ) {
            GetDocOutput()->AdjustSourcePath(fLineComment, relpath);
            out << fLineComment << std::endl;
         } else if (!InContext(kDirective)) {
            GetDocOutput()->AdjustSourcePath(fLineSource, relpath);
            out << fLineSource << std::endl;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand keywords in text, writing to out.

void TDocParser::DecorateKeywords(std::ostream& out, const char *text)
{
   TString str(text);
   DecorateKeywords(str);
   out << str;
}

////////////////////////////////////////////////////////////////////////////////
/// Find keywords in line and create URLs around them. Escape characters with a
/// special meaning for HTML. Protect "Begin_Html"/"End_Html" pairs, and set the
/// parsing context. Evaluate sequences like a::b->c.
/// Skip regions where directives are active.

void TDocParser::DecorateKeywords(TString& line)
{
   std::list<TClass*> currentType;

   enum {
      kNada,
      kMember,
      kScope,
      kNumAccesses
   } scoping = kNada;

   currentType.push_back(0);

   Ssiz_t i = 0;
   while (isspace((UChar_t)line[i]))
      ++i;

   Ssiz_t startOfLine = i;

   // changed when the end of a directive is encountered, i.e.
   // from where fLineSource needs to be appended to fLineComment
   Ssiz_t copiedToCommentUpTo = 0;

   if (InContext(kDirective) && fDirectiveHandlers.Last()) {
      // we're only waiting for an "End_Whatever" and ignoring everything else
      TDocDirective* directive = (TDocDirective*)fDirectiveHandlers.Last();
      const char* endTag = directive->GetEndTag();
      Ssiz_t posEndTag = i;
      while (kNPOS != (posEndTag = line.Index(endTag, posEndTag, TString::kIgnoreCase)))
         if (posEndTag == 0 || line[posEndTag - 1] != '"') // escaping '"'
            break;
      if (posEndTag != kNPOS)
         i = posEndTag;
      else {
         Ssiz_t start = 0;
         if (!InContext(kComment) || (InContext(kComment) & kCXXComment)) {
            // means we are in a C++ comment
            while (isspace((UChar_t)fLineRaw[start])) ++start;
            if (fLineRaw[start] == '/' && fLineRaw[start + 1] == '/')
               start += 2;
            else start = 0;
         }
         directive->AddLine(fLineRaw(start, fLineRaw.Length()));
         while(i < line.Length())
            fDocOutput->ReplaceSpecialChars(line, i);
         copiedToCommentUpTo = i;
      }
   }

   for (; i < line.Length(); ++i) {

      if (!currentType.back())
         scoping = kNada;

      // evaluate scope relation
      if (Context() == kCode
         || Context() == kComment) {
         if (currentType.back())
            switch (line[i]) {
               case ':':
                  if (line[i + 1] == ':') {
                     scoping = kScope;
                     i += 1;
                     continue;
                  }
                  break;
               case '-':
                  if (line[i + 1] == '>') {
                     scoping = kMember;
                     i += 1;
                     continue;
                  }
                  break;
               case '.':
                  if (line[i + 1] != '.') {
                     // prevent "..."
                     scoping = kMember;
                     continue;
                  }
                  break;
            }
         switch (line[i]) {
            case '(':
               currentType.push_back(0);
               scoping = kNada;
               continue;
               break;
            case ')':
               if (currentType.size() > 1)
                  currentType.pop_back();
               scoping = kMember;
               continue;
               break;
         }
         if (i >= line.Length())
            break;
      } else // code or comment
         currentType.back() = 0;


      if (!IsWord(line[i])){

         Bool_t haveHtmlEscapedChar = Context() == kString
            && i > 2 && line[i] == '\'' && line[i-1] == ';';
         if (haveHtmlEscapedChar) {
            Ssiz_t posBegin = i - 2;
            while (posBegin > 0 && IsWord(line[posBegin]))
               --posBegin;
            haveHtmlEscapedChar = posBegin > 0 &&
               line[posBegin] == '&' && line[posBegin - 1] == '\'';
         }
         EParseContext context = Context();
         Bool_t closeString = context == kString
            && (  line[i] == '"'
               || (line[i] == '\''
                   && (  (i > 1 && line[i - 2] == '\'')
                      || (i > 3 && line[i - 2] == '\\' && line[i - 3] == '\'')))
               || haveHtmlEscapedChar)
            && (i == 0 || line[i - 1] != '\\'); // but not "foo \"str...
         if (context == kCode || context == kComment) {
            if (line[i] == '"' || (line[i] == '\'' && (
                  // 'a'
                  (line.Length() > i + 2 && line[i + 2] == '\'') ||
                  // '\a'
                  (line.Length() > i + 3 && line[i + 1] == '\'' && line[i + 3] == '\'')))) {

               fDocOutput->DecorateEntityBegin(line, i, kString);
               fParseContext.push_back(kString);
               currentType.back() = 0;
               closeString = kFALSE;
            } else if (context == kCode
               && line[i] == '/' && (line[i+1] == '/' || line[i+1] == '*')) {
               fParseContext.push_back(kComment);
               if (line[i+1] == '/')
                  fParseContext.back() |= kCXXComment;
               currentType.back() = 0;
               fDocOutput->DecorateEntityBegin(line, i, kComment);
               ++i;
            } else if (context == kComment
               && !(fParseContext.back() & kCXXComment)
               && line.Length() > i + 1
               && line[i] == '*' && line[i+1] == '/') {
               if (fParseContext.size()>1)
                  fParseContext.pop_back();

               currentType.back() = 0;
               i += 2;
               fDocOutput->DecorateEntityEnd(line, i, kComment);
               if (!fCommentAtBOL) {
                  if (InContext(kDirective))
                     ((TDocDirective*)fDirectiveHandlers.Last())->AddLine(line(copiedToCommentUpTo, i));
                  else
                     fLineComment += line(copiedToCommentUpTo, i);
                  copiedToCommentUpTo = i;
               }
            } else if (startOfLine == i
               && line[i] == '#'
               && context == kCode) {
               ExpandCPPLine(line, i);
            }
         } // if context is comment or code

         if (i < line.Length())
            fDocOutput->ReplaceSpecialChars(line, i);

         if (closeString) {
            fDocOutput->DecorateEntityEnd(line, i, kString);
            if (fParseContext.size()>1)
               fParseContext.pop_back();

            currentType.back() = 0;
         }
         --i; // i already moved by ReplaceSpecialChar

         continue;
      } // end of "not a word"

      // get the word
      Ssiz_t endWord = i;
      while (endWord < line.Length() && IsName(line[endWord]))
         endWord++;

      if (Context() == kString || Context() == kCPP) {
         // don't replace in strings, cpp, etc
         i = endWord - 1;
         continue;
      }

      TString word(line(i, endWord - i));

      // '"' escapes handling of "Begin_..."/"End_..."
      if ((i == 0 || (i > 0 && line[i - 1] != '"'))
         && HandleDirective(line, i, word, copiedToCommentUpTo)) {
         // something special happened; the currentType is gone.
         currentType.back() = 0;
         continue;
      }

      // don't replace keywords in comments
      if (Context() == kCode
         && fgKeywords.find(word.Data()) != fgKeywords.end()) {
         fDocOutput->DecorateEntityBegin(line, i, kKeyword);
         i += word.Length();
         fDocOutput->DecorateEntityEnd(line, i, kKeyword);
         --i; // -1 for ++i
         currentType.back() = 0;
         continue;
      }

      // Now decorate scopes and member, referencing their documentation:

      // generic layout:
      // A::B::C::member[arr]->othermember
      // we iterate through this, first scope is A, and currentType will be set toA,
      // next we see ::B, "::" signals to use currentType,...

      TDataType* subType = 0;
      TClass* subClass = 0;
      TDataMember *datamem = 0;
      TMethod *meth = 0;
      const char* globalTypeName = 0;
      if (currentType.empty()) {
         Warning("DecorateKeywords", "type context is empty!");
         currentType.push_back(0);
      }
      TClass* lookupScope = currentType.back();

      if (scoping == kNada) {
         if (fCurrentClass)
            lookupScope = fCurrentClass;
         else
            lookupScope = fRecentClass;
      }

      if (scoping == kNada) {
         subType = gROOT->GetType(word);
         if (!subType)
            subClass = fHtml->GetClass(word);
         if (!subType && !subClass) {
            TGlobal *global = gROOT->GetGlobal(word);
            if (global) {
               // cannot doc globals; take at least their type...
               globalTypeName = global->GetTypeName();
               subClass = fHtml->GetClass(globalTypeName);
               if (!subClass)
                  subType = gROOT->GetType(globalTypeName);
               else // hack to prevent current THtml obj from showing up - we only want gHtml
                  if (subClass == THtml::Class() && word != "gHtml")
                     subClass = 0;
            }
         }
         if (!subType && !subClass) {
            // too bad - cannot doc yet...
            //TFunction *globFunc = gROOT->GetGlobalFunctionWithPrototype(word);
            //globFunc = 0;
         }
         if (!subType && !subClass) {
            // also try template
            while (isspace(line[endWord])) ++endWord;
            if (line[endWord] == '<' || line[endWord] == '>') {
               // check for possible template
               Ssiz_t endWordT = endWord + 1;
               int templateLevel = 1;
               while (endWordT < line.Length()
                      && (templateLevel
                          || IsName(line[endWordT])
                          || line[endWordT] == '<'
                          || line[endWordT] == '>')) {
                  if (line[endWordT] == '<')
                     ++templateLevel;
                  else if (line[endWordT] == '>')
                     --templateLevel;
                  endWordT++;
               }
               subClass = fHtml->GetClass(line(i, endWordT - i).Data());
               if (subClass)
                  word = line(i, endWordT - i);
            }
         }
      }

      if (lookupScope && !subType && !subClass) {
         if (scoping == kScope) {
            TString subClassName(lookupScope->GetName());
            subClassName += "::";
            subClassName += word;
            subClass = fHtml->GetClass(subClassName);
            if (!subClass)
               subType = gROOT->GetType(subClassName);
         }
         if (!subClass && !subType) {
            // also try A::B::c()
            datamem = lookupScope->GetDataMember(word);
            if (!datamem)
               meth = lookupScope->GetMethodAllAny(word);
         }
         if (!subClass && !subType && !datamem && !meth) {
            // also try template
            while (isspace(line[endWord])) ++endWord;
            if (line[endWord] == '<' || line[endWord] == '>') {
               // check for possible template
               Ssiz_t endWordT = endWord + 1;
               int templateLevel = 1;
               while (endWordT < line.Length()
                      && (templateLevel
                          || IsName(line[endWordT])
                          || line[endWordT] == '<'
                          || line[endWordT] == '>')) {
                  if (line[endWordT] == '<')
                     ++templateLevel;
                  else if (line[endWordT] == '>')
                     --templateLevel;
                  endWordT++;
               }
               TString subClassName(lookupScope->GetName());
               subClassName += "::";
               subClassName += line(i, endWordT - i);
               subClass = fHtml->GetClass(subClassName);
               if (subClass)
                  word = line(i, endWordT - i);
            }
         }
      }
      // create the link
      TString mangledWord(word);
      fDocOutput->ReplaceSpecialChars(mangledWord);
      line.Replace(i, word.Length(), mangledWord);

      TSubString substr(line(i, mangledWord.Length()));
      if (subType) {
         fDocOutput->ReferenceEntity(substr, subType,
            globalTypeName ? globalTypeName : subType->GetName());
         currentType.back() = 0;
      } else if (subClass) {
         fDocOutput->ReferenceEntity(substr, subClass,
            globalTypeName ? globalTypeName : subClass->GetName());

         currentType.back() = subClass;
         fRecentClass = subClass;
      } else if (datamem || meth) {
            if (datamem) {
               fDocOutput->ReferenceEntity(substr, datamem);

               if (datamem->GetTypeName())
                  currentType.back() = fHtml->GetClass(datamem->GetTypeName());
            } else {
               fDocOutput->ReferenceEntity(substr, meth);

               TString retTypeName = meth->GetReturnTypeName();
               if (retTypeName.BeginsWith("const "))
                  retTypeName.Remove(0,6);
               Ssiz_t pos=0;
               while (IsWord(retTypeName[pos]) || retTypeName[pos]=='<' || retTypeName[pos]=='>' || retTypeName[pos]==':')
                  ++pos;
               retTypeName.Remove(pos, retTypeName.Length());
               if (retTypeName.Length())
                  currentType.back() = fHtml->GetClass(retTypeName);
            }
      } else
         currentType.back() = 0;

      //i += mangledWord.Length();
      i += substr.Length();

      --i; // due to ++i
   } // while i < line.Length()
   if (i > line.Length())
      i = line.Length();

   // clean up, no strings across lines
   if (Context() == kString) {
      fDocOutput->DecorateEntityEnd(line, i, kString);
      if (fParseContext.size()>1)
         fParseContext.pop_back();
      currentType.back() = 0;
   }

   // HandleDirective already copied the chunk before the directive
   // from fLineSource to fLineComment. So we're done up to "i" in
   // fLineSource; next time we encounter a directive we just need
   // to copy from startOfComment on.
   if ((InContext(kComment) || fCommentAtBOL) && copiedToCommentUpTo < line.Length()) {
      if (InContext(kDirective))
         ((TDocDirective*)fDirectiveHandlers.Last())->AddLine(line(copiedToCommentUpTo, line.Length()));
      else
         fLineComment += line(copiedToCommentUpTo, line.Length());
   }

   // Do this after we append to fLineComment, otherwise the closing
   // </span> gets sent to the directive.
   // clean up, no CPP comment across lines
   if (InContext(kComment) & kCXXComment) {
      fDocOutput->DecorateEntityEnd(line, i, kComment);
      if (fLineComment.Length()) {
         Ssiz_t pos = fLineComment.Length();
         fDocOutput->DecorateEntityEnd(fLineComment, pos, kComment);
      }
      RemoveCommentContext(kTRUE);
      currentType.back() = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// reduce method count for method called name,
/// removing it from fMethodCounts once the count reaches 0.

void TDocParser::DecrementMethodCount(const char* name)
{
   typedef std::map<std::string /*method name*/, Int_t > MethodCount_t;
   MethodCount_t::iterator iMethodName = fMethodCounts.find(name);
   if (iMethodName != fMethodCounts.end()) {
      --(iMethodName->second);
      if (iMethodName->second <= 0)
         fMethodCounts.erase(iMethodName);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete output generated by prior runs of all known directives;
/// the output file names might have changes.

void  TDocParser::DeleteDirectiveOutput() const
{
   TIter iClass(gROOT->GetListOfClasses());
   TClass* cl = 0;
   while ((cl = (TClass*) iClass()))
      if (cl != TDocDirective::Class()
         && cl->InheritsFrom(TDocDirective::Class())) {
         TDocDirective* directive = (TDocDirective*) cl->New();
         if (!directive) continue;
         directive->SetParser(const_cast<TDocParser*>(this));
         directive->DeleteOutput();
         delete directive;
      }
}

////////////////////////////////////////////////////////////////////////////////
/// Expand preprocessor statements
///
///
/// Input: line - line containing the CPP statement,
///        pos  - position of '#'
///
///  NOTE: Looks for the #include statements and
///        creates link to the corresponding file
///        if such file exists
///

void TDocParser::ExpandCPPLine(TString& line, Ssiz_t& pos)
{
   Bool_t linkExist    = kFALSE;
   Ssiz_t posEndOfLine = line.Length();
   Ssiz_t posHash      = pos;

   Ssiz_t posInclude = line.Index("include", pos);
   if (posInclude != kNPOS) {
      TString filename;
      Ssiz_t posStartFilename = posInclude + 7;
      if (line.Tokenize(filename, posStartFilename, "[<\"]")) {
         Ssiz_t posEndFilename = posStartFilename;
         if (line.Tokenize(filename, posEndFilename, "[>\"]")) {
            R__LOCKGUARD(fHtml->GetMakeClassMutex());

            TString filesysFileName;
            if (fHtml->GetPathDefinition().GetFileNameFromInclude(filename, filesysFileName)) {
               fDocOutput->CopyHtmlFile(filesysFileName);

               TString endOfLine(line(posEndFilename - 1, line.Length()));
               line.Remove(posStartFilename, line.Length());
               for (Ssiz_t i = pos; i < line.Length();)
                  fDocOutput->ReplaceSpecialChars(line, i);

               line += "<a href=\"./";
               line += gSystem->BaseName(filename);
               line += "\">";
               line += filename + "</a>" + endOfLine[0]; // add include file's closing '>' or '"'
               posEndOfLine = line.Length() - 1; // set the "processed up to" to it
               fDocOutput->ReplaceSpecialChars(line, posEndOfLine); // and run replace-special-char on it

               line += endOfLine(1, endOfLine.Length()); // add the unprocessed part of the line back

               linkExist = kTRUE;
            }
         }
      }
   }

   if (!linkExist) {
      fDocOutput->ReplaceSpecialChars(line);
      posEndOfLine = line.Length();
   }

   Ssiz_t posHashAfterDecoration = posHash;
   fDocOutput->DecorateEntityBegin(line, posHashAfterDecoration, kCPP);
   posEndOfLine += posHashAfterDecoration - posHash;

   fDocOutput->DecorateEntityEnd(line, posEndOfLine, kCPP);
   pos = posEndOfLine;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the name of module for which sources are currently parsed.

void TDocParser::GetCurrentModule(TString& out_module) const {
   if (fCurrentModule) out_module = fCurrentModule;
   else if (fCurrentClass) fHtml->GetModuleNameForClass(out_module, fCurrentClass);
   else out_module = "(UNKNOWN MODULE WHILE PARSING)";
}

////////////////////////////////////////////////////////////////////////////////
/// Process directives to the documentation engine, like "Begin_Html" / "End_Html",
/// "Begin_Macro" / "End_Macro", and "Begin_Latex" / "End_Latex".

Bool_t TDocParser::HandleDirective(TString& line, Ssiz_t& pos, TString& word,
                                   Ssiz_t& copiedToCommentUpTo)
{
   Bool_t begin = kTRUE;
   TClass* clDirective = IsDirective(line, pos, word, begin);
   if (!clDirective)
      return kFALSE;

   // we'll need end later on: afer the begin block, both end _and_ begin can be true.
   Bool_t end = !begin;

   TDocDirective* directive = 0; // allow re-use of object from begin block in end

   if (begin) {
      // copy from fLineSource to fLineComment, starting at copiedToCommentUpTo
      if (InContext(kDirective))
         ((TDocDirective*)fDirectiveHandlers.Last())->AddLine(fLineSource(copiedToCommentUpTo, pos - copiedToCommentUpTo));
      else
         fLineComment += fLineSource(copiedToCommentUpTo, pos - copiedToCommentUpTo);
      copiedToCommentUpTo = pos;

      pos += word.Length(); // skip the keyword

      directive = (TDocDirective*) clDirective->New();
      if (!directive)
         return kFALSE;

      directive->SetParser(this);
      if (fCurrentMethodTag.Length())
         directive->SetTag(fCurrentMethodTag);
      directive->SetCounter(fDirectiveCount++);

      // parse parameters
      TString params;
      if (begin && line[pos] == '(') {
         std::list<char> waitForClosing;
         Ssiz_t endParam = pos + 1;
         for (; endParam < line.Length()
            && (line[endParam] != ')' || !waitForClosing.empty()); ++endParam) {
            const char c = line[endParam];
            if (!waitForClosing.empty() && waitForClosing.back() == c) {
               waitForClosing.pop_back();
               continue;
            }
            switch (c) {
               case '"':
                  if (waitForClosing.empty() || waitForClosing.back() != '\'')
                     waitForClosing.push_back('"');
                  break;
               case '\'':
                  if (waitForClosing.empty() || waitForClosing.back() != '"')
                     waitForClosing.push_back('\'');
                  break;
               case '(':
                  if (waitForClosing.empty() || (waitForClosing.back() != '"' && waitForClosing.back() != '\''))
                     waitForClosing.push_back(')');
                  break;
               case '\\':
                  ++endParam; // skip next char
               default:
                  break;
            };
         }
         if (waitForClosing.empty()) {
            params = line(pos + 1, endParam - (pos + 1));
            pos += params.Length() + 2; // params + parentheses
         }
         directive->SetParameters(params);
      }

      // check for end tag in current line
      Ssiz_t posEndTag = pos;
      const char* endTag = directive->GetEndTag();
      Ssiz_t lenEndTag = strlen(endTag);
      while (kNPOS != (posEndTag = line.Index(endTag, posEndTag, TString::kIgnoreCase))) {
         if (line[posEndTag - 1] == '"') {
            posEndTag += lenEndTag;
            continue; // escaping '"'
         }
         break;
      }
      if (posEndTag != kNPOS) {
         end = kTRUE; // we just continue below!
      } else {
         fDirectiveHandlers.AddLast(directive);

         fParseContext.push_back(kDirective);
         if (InContext(kComment) & kCXXComment)
            fParseContext.back() |= kCXXComment;

         posEndTag = line.Length();
      }

      directive->AddLine(line(pos, posEndTag - pos));
      TString remainder(line(posEndTag, line.Length()));
      line.Remove(posEndTag, line.Length());

      while (pos < line.Length())
         fDocOutput->ReplaceSpecialChars(line, pos);

      pos = line.Length();
      // skip the remainder of the line
      copiedToCommentUpTo = line.Length();
      line += remainder;
   }

   // no else - "end" can also be set by begin having an end tag!
   if (end) {

      if (!begin)
         pos += word.Length(); // skip the keyword
      else pos += word.Length() - 2; // "Begin" is 2 chars longer than "End"

      if (!directive) directive = (TDocDirective*) fDirectiveHandlers.Last();

      if (!directive) {
         Warning("HandleDirective", "Cannot find directive handler object %s !",
            fLineRaw.Data());
         return kFALSE;
      }

      if (!begin) {
         Ssiz_t start = 0;
         if (!InContext(kComment) || (InContext(kComment) & kCXXComment)) {
            // means we are in a C++ comment
            while (isspace((UChar_t)fLineRaw[start])) ++start;
            if (fLineRaw[start] == '/' && fLineRaw[start + 1] == '/')
               start += 2;
            else start = 0;
         }
         directive->AddLine(line(start, pos - word.Length() - start));

         TString remainder(line(pos, line.Length()));
         line.Remove(pos, line.Length());
         fDocOutput->ReplaceSpecialChars(line);
         pos = line.Length();
         line += remainder;
      }
      copiedToCommentUpTo = pos;

      TString result;
      directive->GetResult(result);

      if (!begin)
         fDirectiveHandlers.Remove(fDirectiveHandlers.LastLink());
      delete directive;

      if (!begin) {
         // common to all directives: pop context
         Bool_t isInCxxComment = InContext(kDirective) & kCXXComment;
         if (fParseContext.size()>1)
            fParseContext.pop_back();
         if (isInCxxComment && !InContext(kComment)) {
            fParseContext.push_back(kComment | kCXXComment);
            fDocOutput->DecorateEntityBegin(line, pos, kComment);
         }
      }

      if (InContext(kDirective) && fDirectiveHandlers.Last())
         ((TDocDirective*)fDirectiveHandlers.Last())->AddLine(result(0, result.Length()));
      else
         fLineComment += result;

      /* NO - this can happen e.g. for "BEGIN_HTML / *..." (see doc in this class)
      if (Context() != kComment) {
         Warning("HandleDirective", "Popping back a directive context, but enclosing context is not a comment! At:\n%s",
            fLineRaw.Data());
         fParseContext.push_back(kComment);
      }
      */
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// checks whether we are in a parse context, return the entry closest
/// to the current context.
/// If context is a EParseContextFlag just look for the first match in
/// the flags

UInt_t TDocParser::InContext(Int_t context) const
{
   UInt_t lowerContext = context & kParseContextMask;
   UInt_t contextFlag  = context & kParseContextFlagMask;

   for (std::list<UInt_t>::const_reverse_iterator iPC = fParseContext.rbegin();
      iPC != fParseContext.rend(); ++iPC)
      if (!lowerContext || ((lowerContext && ((*iPC & kParseContextMask) == lowerContext))
         && (!contextFlag || (contextFlag && (*iPC & contextFlag)))))
         return *iPC;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// fill C++ keywords into fgKeywords

void TDocParser::InitKeywords() const
{
   if (!fgKeywords.empty())
      return;

   fgKeywords.insert("asm");
   fgKeywords.insert("auto");
   fgKeywords.insert("bool");
   fgKeywords.insert("break");
   fgKeywords.insert("case");
   fgKeywords.insert("catch");
   fgKeywords.insert("char");
   fgKeywords.insert("class");
   fgKeywords.insert("const");
   fgKeywords.insert("const_cast");
   fgKeywords.insert("continue");
   fgKeywords.insert("default");
   fgKeywords.insert("delete");
   fgKeywords.insert("do");
   fgKeywords.insert("double");
   fgKeywords.insert("dynamic_cast");
   fgKeywords.insert("else");
   fgKeywords.insert("enum");
   fgKeywords.insert("explicit");
   fgKeywords.insert("export");
   fgKeywords.insert("extern");
   fgKeywords.insert("false");
   fgKeywords.insert("float");
   fgKeywords.insert("for");
   fgKeywords.insert("friend");
   fgKeywords.insert("goto");
   fgKeywords.insert("if");
   fgKeywords.insert("inline");
   fgKeywords.insert("int");
   fgKeywords.insert("long");
   fgKeywords.insert("mutable");
   fgKeywords.insert("namespace");
   fgKeywords.insert("new");
   fgKeywords.insert("operator");
   fgKeywords.insert("private");
   fgKeywords.insert("protected");
   fgKeywords.insert("public");
   fgKeywords.insert("register");
   fgKeywords.insert("reinterpret_cast");
   fgKeywords.insert("return");
   fgKeywords.insert("short");
   fgKeywords.insert("signed");
   fgKeywords.insert("sizeof");
   fgKeywords.insert("static");
   fgKeywords.insert("static_cast");
   fgKeywords.insert("struct");
   fgKeywords.insert("switch");
   fgKeywords.insert("template");
   fgKeywords.insert("this");
   fgKeywords.insert("throw");
   fgKeywords.insert("true");
   fgKeywords.insert("try");
   fgKeywords.insert("typedef");
   fgKeywords.insert("typeid");
   fgKeywords.insert("typename");
   fgKeywords.insert("union");
   fgKeywords.insert("unsigned");
   fgKeywords.insert("using");
   fgKeywords.insert("virtual");
   fgKeywords.insert("void");
   fgKeywords.insert("volatile");
   fgKeywords.insert("wchar_t");
   fgKeywords.insert("while");
}

////////////////////////////////////////////////////////////////////////////////
/// return whether word at line's pos is a valid directive, and returns its
/// TDocDirective's TClass object, or 0 if it's not a directive. Set begin
/// to kTRUE for "Begin_..."
/// You can implement your own handlers by implementing a class deriving
/// from TDocHandler, and calling it TDocTagDirective for "BEGIN_TAG",
/// "END_TAG" blocks.

TClass* TDocParser::IsDirective(const TString& line, Ssiz_t pos,
                                        const TString& word, Bool_t& begin) const
{
   // '"' serves as escape char
   if (pos > 0 &&  line[pos - 1] == '"')
      return 0;

   begin      = word.BeginsWith("begin_", TString::kIgnoreCase);
   Bool_t end = word.BeginsWith("end_", TString::kIgnoreCase);

   if (!begin && !end)
      return 0;

   /* NO - we can have "BEGIN_HTML / * ..."
   if (!InContext(kComment))
      return 0;
   */

   TString tag = word( begin ? 6 : 4, word.Length());

   if (!tag.Length())
      return 0;

   tag.ToLower();
   tag[0] -= 'a' - 'A'; // first char is caps
   tag.Prepend("TDoc");
   tag += "Directive";

   TClass* clDirective = TClass::GetClass(tag, kFALSE);

   if (gDebug > 0 && !clDirective)
      Warning("IsDirective", "Unknown THtml directive %s in line %d!", word.Data(), fLineNo);

   return clDirective;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if c is a valid C++ name character
///
///
///  Input: c - a single character
///
/// Output: TRUE if c is a valid C++ name character
///         and FALSE if it's not.
///
///   NOTE: Valid name characters are [a..zA..Z0..9_~],
///

Bool_t TDocParser::IsName(UChar_t c)
{
   Bool_t ret = kFALSE;

   if (isalnum(c) || c == '_' || c == '~')
      ret = kTRUE;

   return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// Check if c is a valid first character for C++ name
///
///
///  Input: c - a single character
///
/// Output: TRUE if c is a valid first character for C++ name,
///         and FALSE if it's not.
///
///   NOTE: Valid first characters are [a..zA..Z_~]
///

Bool_t TDocParser::IsWord(UChar_t c)
{
   Bool_t ret = kFALSE;

   if (isalpha(c) || c == '_' || c == '~')
      ret = kTRUE;

   return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// Search for a method starting at posMethodName, and return its return type,
/// its name, and its arguments. If the end of arguments is not found in the
/// current line, get a new line from sourceFile, beautify it to srcOut, creating
/// an anchor as necessary. When this function returns, posMethodName points to the
/// end of the function declaration, i.e. right after the arguments' closing bracket.
/// If posMethodName == kNPOS, we look for the first matching method in fMethodCounts.

TMethod* TDocParser::LocateMethodInCurrentLine(Ssiz_t &posMethodName, TString& ret,
                                               TString& name, TString& params,
                                               Bool_t& isconst, std::ostream &srcOut,
                                               TString &anchor, std::ifstream& sourceFile,
                                               Bool_t allowPureVirtual)
{
   typedef std::map<std::string /*method name*/, Int_t > MethodCount_t;
   isconst = false;

   if (posMethodName == kNPOS) {
      name.Remove(0);
      TMethod * meth = 0;
      Ssiz_t posBlock = fLineRaw.Index('{');
      Ssiz_t posQuote = fLineRaw.Index('"');
      if (posQuote != kNPOS && (posBlock == kNPOS || posQuote < posBlock))
         posBlock = posQuote;
      if (posBlock == kNPOS)
         posBlock = fLineRaw.Length();
      for (MethodCount_t::iterator iMethodName = fMethodCounts.begin();
         !name.Length() && iMethodName != fMethodCounts.end(); ++iMethodName) {
         TString lookFor(iMethodName->first);
         posMethodName = fLineRaw.Index(lookFor);
         if (posMethodName != kNPOS && posMethodName < posBlock
            && (posMethodName == 0 || !IsWord(fLineRaw[posMethodName - 1]))) {
            // check whether the method name is followed by optional spaces and
            // an opening parathesis
            Ssiz_t posMethodEnd = posMethodName + lookFor.Length();
            while (isspace((UChar_t)fLineRaw[posMethodEnd])) ++posMethodEnd;
            if (fLineRaw[posMethodEnd] == '(') {
               meth = LocateMethodInCurrentLine(posMethodName, ret, name, params, isconst,
                                                srcOut, anchor, sourceFile, allowPureVirtual);
               if (name.Length())
                  return meth;
            }
         }
      }
      return 0;
   }

   name = fLineRaw(posMethodName, fLineRaw.Length() - posMethodName);

   // extract return type
   ret = fLineRaw(0, posMethodName);
   if (ret.Length()) {
      while (ret.Length() && (IsName(ret[ret.Length() - 1]) || ret[ret.Length()-1] == ':'))
         ret.Remove(ret.Length() - 1, 1);
      Strip(ret);
      Bool_t didSomething = kTRUE;
      while (didSomething) {
         didSomething = kFALSE;
         if (ret.BeginsWith("inline ")) {
            didSomething = kTRUE;
            ret.Remove(0, 7);
         }
         if (ret.BeginsWith("static ")) {
            didSomething = kTRUE;
            ret.Remove(0, 7);
         }
         if (ret.BeginsWith("virtual ")) {
            didSomething = kTRUE;
            ret.Remove(0, 8);
         }
      } // while replacing static, virtual, inline
      Strip(ret);
   }

   // extract parameters
   Ssiz_t posParam = name.First('(');
   if (posParam == kNPOS ||
      // no strange return types, please
      ret.Contains("{") || ret.Contains("}") || ret.Contains("(") || ret.Contains(")")
      || ret.Contains("=")) {
      ret.Remove(0);
      name.Remove(0);
      params.Remove(0);
      return 0;
   }

   if (name.BeginsWith("operator")) {
      // op () (...)
      Ssiz_t checkOpBracketParam = posParam + 1;
      while (isspace((UChar_t)name[checkOpBracketParam]))
         ++checkOpBracketParam;
      if (name[checkOpBracketParam] == ')') {
         ++checkOpBracketParam;
         while (isspace((UChar_t)name[checkOpBracketParam]))
            ++checkOpBracketParam;
         if (name[checkOpBracketParam] == '(')
            posParam = checkOpBracketParam;
      }
   } // check for op () (...)

   if (posParam == kNPOS) {
      ret.Remove(0);
      name.Remove(0);
      params.Remove(0);
      return 0;
   }

   params = name(posParam, name.Length() - posParam);
   name.Remove(posParam);
   while (name.Length() && isspace((UChar_t)name[name.Length() - 1]))
      name.Remove(name.Length() - 1);
   if (!name.Length()) {
      ret.Remove(0);
      name.Remove(0);
      params.Remove(0);
      return 0;
   }

   MethodCount_t::const_iterator iMethodName = fMethodCounts.find(name.Data());
   if (iMethodName == fMethodCounts.end() || iMethodName->second <= 0) {
      ret.Remove(0);
      name.Remove(0);
      params.Remove(0);
      return 0;
   }

   // find end of param
   Ssiz_t posParamEnd = 1;
   Int_t bracketLevel = 1;
   while (bracketLevel) {
      const char* paramEnd = strpbrk(params.Data() + posParamEnd, ")(\"'");
      if (!paramEnd) {
         // func with params over multiple lines
         // gotta write out this line before it gets lost
         if (!anchor.Length()) {
            // request an anchor, just in case...
            AnchorFromLine(fLineStripped, anchor);
            if (srcOut)
               srcOut << "<a name=\"" << anchor << "\"></a>";
         }
         ++fLineNumber;
         if (srcOut)
            WriteSourceLine(srcOut);

         fLineRaw.ReadLine(sourceFile, kFALSE);
         if (sourceFile.eof()) {
            Error("LocateMethodInCurrentLine",
               "Cannot find end of signature for function %s!",
               name.Data());
            break;
         }

         fCommentAtBOL = kFALSE;

         // replace class names etc
         fLineStripped = fLineRaw;
         Strip(fLineStripped);

         fLineSource = fLineRaw;
         DecorateKeywords(fLineSource);

         posParamEnd = params.Length();
         params += fLineRaw;
      } else
         posParamEnd = paramEnd - params.Data();
      switch (params[posParamEnd]) {
         case '(': ++bracketLevel; ++posParamEnd; break;
         case ')': --bracketLevel; ++posParamEnd; break;
         case '"': // skip ")"
            ++posParamEnd;
            while (params.Length() > posParamEnd && params[posParamEnd] != '"') {
               // skip '\"'
               if (params[posParamEnd] == '\\') ++posParamEnd;
               ++posParamEnd;
            }
            if (params.Length() <= posParamEnd) {
               // something is seriously wrong - skip :-/
               ret.Remove(0);
               name.Remove(0);
               params.Remove(0);
               return 0;
            }
            ++posParamEnd; // skip trailing '"'
            break;
         case '\'': // skip ')'
            ++posParamEnd;
            if (params[posParamEnd] == '\\') ++posParamEnd;
            posParamEnd += 2;
            break;
         default:
            ++posParamEnd;
      }
   } // while bracketlevel, i.e. (...(..)...)

   {
      TString pastParams(params(posParamEnd, params.Length()));
      pastParams = pastParams.Strip(TString::kLeading);
      isconst = pastParams.BeginsWith("const") && !(isalnum(pastParams[5]) || pastParams[5] == '_');
   }

   Ssiz_t posBlock     = params.Index('{', posParamEnd);
   Ssiz_t posSemicolon = params.Index(';', posParamEnd);
   Ssiz_t posPureVirt  = params.Index('=', posParamEnd);
   if (posSemicolon != kNPOS)
      if ((posBlock == kNPOS || (posSemicolon < posBlock)) &&
         (posPureVirt == kNPOS || !allowPureVirtual)
         && !allowPureVirtual) // allow any "func();" if pv is allowed
         params.Remove(0);

   if (params.Length())
      params.Remove(posParamEnd);

   if (!params.Length()) {
      ret.Remove(0);
      name.Remove(0);
      return 0;
   }
   // update posMethodName to point behind the method
   posMethodName = posParam + posParamEnd;
   if (fCurrentClass) {
      TMethod* meth = fCurrentClass->GetMethodAny(name);
      if (meth) {
         fDirectiveCount = 0;
         fCurrentMethodTag = name + "_";
         fCurrentMethodTag += fMethodCounts[name.Data()];
         return meth;
      }
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Locate methods, starting in the source file, then inline, then
/// immediately inside the class declaration. While doing that also
/// find the class description and special tags like the macro tag etc.

void TDocParser::Parse(std::ostream& out)
{
   fClassDocState = kClassDoc_LookingNothingFound;

   DeleteDirectiveOutput();

   LocateMethodsInSource(out);
   LocateMethodsInHeaderInline(out);
   LocateMethodsInHeaderClassDecl(out);

   if (!fSourceInfo[kInfoLastUpdate].Length()) {
      TDatime date;
      fSourceInfo[kInfoLastUpdate] = date.AsString();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Collect methods from the source or header file called filename.
/// It generates a beautified version of the source file on the fly;
/// the output file is given by the fCurrentClass's name, and sourceExt.
/// Documentation is extracted to out.
///   lookForSourceInfo: if set, author, lastUpdate, and copyright are
///     extracted (i.e. the values contained in fSourceInfo)
///   useDocxxStyle: if set, documentation can be in front of the method
///     name, not only inside the method. Useful doc Doc++/Doxygen style,
///     and inline methods.
///   lookForClassDescr: if set, the first line matching the class description
///     rules is assumed to be the class description for fCurrentClass; the
///     description is written to out.
///   methodPattern: if set, methods have to be prepended by this tag. Usually
///     the class name + "::". In header files, looking for in-place function
///     definitions, this should be 0. In that case, only functions in
///     fMethodCounts are searched for.

void TDocParser::LocateMethods(std::ostream& out, const char* filename,
                          Bool_t lookForSourceInfo /*= kTRUE*/,
                          Bool_t useDocxxStyle /*= kFALSE*/,
                          Bool_t allowPureVirtual /*= kFALSE*/,
                          const char* methodPattern /*= 0*/,
                          const char* sourceExt /*= 0 */)
{
   TString sourceFileName(filename);
   fCurrentFile = filename;
   if (!sourceFileName.Length()) {
      fHtml->GetImplFileName(fCurrentClass, kFALSE, sourceFileName);
      Error("LocateMethods", "Can't find source file '%s' for class %s!",
         sourceFileName.Data(), fCurrentClass->GetName());
      return;
   }
   std::ifstream sourceFile(sourceFileName.Data());
   if (!sourceFile || !sourceFile.good()) {
      Error("LocateMethods", "Can't open file '%s' for reading!", sourceFileName.Data());
      return;
   }

   TPMERegexp patternRE(methodPattern ? methodPattern : "");

   TString codeOneLiner;
   TString methodRet;
   TString methodName;
   TString methodParam;
   Bool_t methodIsConst = kFALSE;
   TString anchor;
   TString docxxComment;

   Bool_t wroteMethodNowWaitingForOpenBlock = kFALSE;

   std::ofstream srcHtmlOut;
   TString srcHtmlOutName;
   if (sourceExt && sourceExt[0]) {
      static_cast<TClassDocOutput*>(fDocOutput)->CreateSourceOutputStream(srcHtmlOut, sourceExt, srcHtmlOutName);
      fLineNumber = 0;
   } else {
      sourceExt = 0;
      srcHtmlOutName = fCurrentClass->GetName();
      fDocOutput->NameSpace2FileName(srcHtmlOutName);
      gSystem->PrependPathName("src", srcHtmlOutName);
      srcHtmlOutName += ".h.html";
   }

   fParseContext.clear();
   fParseContext.push_back(kCode);
   fDocContext = kIgnore;
   fLineNo = 0;

   while (!sourceFile.eof()) {
      Bool_t needAnchor = kFALSE;

      ++fLineNo; // we count fortrany

      fLineRaw.ReadLine(sourceFile, kFALSE);
      if (sourceFile.eof()) break;

      fCommentAtBOL = InContext(kComment);

      // replace class names etc
      fLineStripped = fLineRaw;
      Strip(fLineStripped);

      fLineSource = fLineRaw;
      fLineComment = "";
      DecorateKeywords(fLineSource);

      if (!ProcessComment()) {
         // not a commented line

         if (fDocContext == kDocClass && fClassDocState < kClassDoc_Written) {
            TString strippedComment(fComment);
            Strip(strippedComment);
            if (strippedComment.Length() > 0) {
               fLastClassDoc = fComment;
               if (fClassDocState == kClassDoc_LookingNothingFound) {
                  fFirstClassDoc = fComment;
                  fClassDocState = kClassDoc_LookingHaveSomething;
               }
            }
            fDocContext = kIgnore;
         }

         Ssiz_t impIdx = fLineStripped.Index("ClassImp(");
         if (impIdx == 0 && fClassDocState == kClassDoc_LookingHaveSomething) {
            TString name(fCurrentClass->GetName());
            // take unscoped version
            Ssiz_t posLastScope = kNPOS;
            while ((posLastScope = name.Index("::")) != kNPOS)
               name.Remove(0, posLastScope + 2);

            Ssiz_t posName = fLineStripped.Index(name, impIdx);
            if (posName != kNPOS) {
               Ssiz_t posClosingParen = posName + name.Length();
               while (isspace(fLineStripped[posClosingParen])) ++posClosingParen;
               if (fLineStripped[posClosingParen] == ')') {
                  WriteClassDoc(out, kFALSE);
                  fDocContext = kIgnore;
               }
            }
         }

         if (fLineStripped.Length())
            // remove last class doc if it not followed by ClassImp
            // (with optional empty lines in between)
            fLastClassDoc = "";

         // write previous method
         if (methodName.Length() && !wroteMethodNowWaitingForOpenBlock) {
            TString savedComment;
            if (useDocxxStyle && docxxComment.Length()) {
               savedComment = fComment;
               fComment = docxxComment;
            }
            WriteMethod(out, methodRet, methodName, methodParam, methodIsConst,
               gSystem->BaseName(srcHtmlOutName), anchor, codeOneLiner);
            docxxComment.Remove(0);
            if (savedComment[0]) {
               fComment = savedComment;
            }
         }

         if (!wroteMethodNowWaitingForOpenBlock) {
            // check for method
            Ssiz_t posPattern = kNPOS;
            if (methodPattern) {
               posPattern = fLineRaw.Index((TPRegexp&)patternRE);
            }
            if (posPattern != kNPOS && methodPattern) {
               // no strings, no blocks in front of function declarations / implementations
               static const char vetoChars[] = "{\"";
               for (int ich = 0; posPattern != kNPOS && vetoChars[ich]; ++ich) {
                  Ssiz_t posVeto = fLineRaw.Index(vetoChars[ich]);
                  if (posVeto != kNPOS && posVeto < posPattern)
                     posPattern = kNPOS;
               }
            }
            if (posPattern != kNPOS || !methodPattern) {
               if (methodPattern) {
                  patternRE.Match(fLineRaw);
                  posPattern += patternRE[0].Length();
               }
               LocateMethodInCurrentLine(posPattern, methodRet, methodName,
                                         methodParam, methodIsConst, srcHtmlOut,
                                         anchor, sourceFile, allowPureVirtual);
               if (methodName.Length()) {
                  fDocContext = kDocFunc;
                  needAnchor = !anchor.Length();
                  if (useDocxxStyle)
                     docxxComment = fComment;
                  fComment.Remove(0);
                  codeOneLiner.Remove(0);

                  wroteMethodNowWaitingForOpenBlock = fLineRaw.Index("{", posPattern) == kNPOS;
                  wroteMethodNowWaitingForOpenBlock &= fLineRaw.Index(";", posPattern) == kNPOS;
               } else if (fLineRaw.First("{};") != kNPOS)
                  // these chars reset the preceding comment
                  fComment.Remove(0);
            } // pattern matches - could be a method
            else
               fComment.Remove(0);
         } else {
            wroteMethodNowWaitingForOpenBlock &= fLineRaw.Index("{") == kNPOS;
            wroteMethodNowWaitingForOpenBlock &= fLineRaw.Index(";") == kNPOS;
         } // if !wroteMethodNowWaitingForOpenBlock

         if (methodName.Length() && !wroteMethodNowWaitingForOpenBlock) {
            // make sure we don't have more '{' in commentLine than in fLineRaw
            if (!codeOneLiner.Length() &&
                fLineSource.CountChar('{') == 1 &&
                fLineSource.CountChar('}') == 1) {
               // a one-liner
               codeOneLiner = fLineSource;
               codeOneLiner.Remove(0, codeOneLiner.Index('{'));
               codeOneLiner.Remove(codeOneLiner.Index('}') + 1);
            }
         } // if method name and '{'
         // else not a comment, and we don't need the previous one:
         else if (!methodName.Length() && !useDocxxStyle)
            fComment.Remove(0);

         if (needAnchor || fExtraLinesWithAnchor.find(fLineNo) != fExtraLinesWithAnchor.end()) {
            AnchorFromLine(fLineStripped, anchor);
            if (sourceExt)
               srcHtmlOut << "<a name=\"" << anchor << "\"></a>";
         }
         // else anchor.Remove(0); - NO! WriteMethod will need it later!
      } // if !comment

      // check for last update,...
      Ssiz_t posTag = kNPOS;
      if (lookForSourceInfo)
         for (Int_t si = 0; si < (Int_t) kNumSourceInfos; ++si)
            if (!fSourceInfo[si].Length() && (posTag = fLineRaw.Index(fSourceInfoTags[si])) != kNPOS) {
               fSourceInfo[si] = fLineRaw(posTag + strlen(fSourceInfoTags[si]), fLineRaw.Length() - posTag);
               if (si == kInfoAuthor)
                  fDocOutput->FixupAuthorSourceInfo(fSourceInfo[kInfoAuthor]);
            }


      // write to .cxx.html
      ++fLineNumber;
      if (srcHtmlOut)
         WriteSourceLine(srcHtmlOut);
      else if (needAnchor)
         fExtraLinesWithAnchor.insert(fLineNo);
   } // while !sourceFile.eof()

   // deal with last func
   if (methodName.Length()) {
      if (useDocxxStyle && docxxComment.Length())
         fComment = docxxComment;
      WriteMethod(out, methodRet, methodName, methodParam, methodIsConst,
         gSystem->BaseName(srcHtmlOutName), anchor, codeOneLiner);
      docxxComment.Remove(0);
   } else
      WriteClassDoc(out);

   srcHtmlOut << "</pre>" << std::endl;

   fDocOutput->WriteLineNumbers(srcHtmlOut, fLineNumber, gSystem->BaseName(fCurrentFile));

   srcHtmlOut << "</div>" << std::endl;

   fDocOutput->WriteHtmlFooter(srcHtmlOut, "../");

   fParseContext.clear();
   fParseContext.push_back(kCode);
   fDocContext = kIgnore;
   fCurrentFile = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Given fCurrentClass, look for methods in its source file,
/// and extract documentation to out, while beautifying the source
/// file in parallel.

void TDocParser::LocateMethodsInSource(std::ostream& out)
{
   // for Doc++ style
   Bool_t useDocxxStyle = (fHtml->GetDocStyle() == "Doc++");

   TString pattern(fCurrentClass->GetName());
   // take unscoped version
   Ssiz_t posLastScope = kNPOS;
   while ((posLastScope = pattern.Index("::")) != kNPOS)
      pattern.Remove(0, posLastScope + 2);
   pattern += "::";

   TString implFileName;
   if (fHtml->GetImplFileName(fCurrentClass, kTRUE, implFileName)) {
      LocateMethods(out, implFileName, kFALSE /*source info*/, useDocxxStyle,
                    kFALSE /*allowPureVirtual*/, pattern, ".cxx.html");
      Ssiz_t posGt = pattern.Index('>');
      if (posGt != kNPOS) {
         // template! Re-run with pattern '...<.*>::'
         Ssiz_t posLt = pattern.Index('<');
         if (posLt != kNPOS && posLt < posGt) {
            pattern.Replace(posLt + 1, posGt - posLt - 1, ".*");
            LocateMethods(out, implFileName, kFALSE /*source info*/, useDocxxStyle,
                    kFALSE /*allowPureVirtual*/, pattern, ".cxx.html");
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Given fCurrentClass, look for methods in its header file,
/// and extract documentation to out.

void TDocParser::LocateMethodsInHeaderInline(std::ostream& out)
{
   // for inline methods, always allow doc before func
   Bool_t useDocxxStyle = kTRUE;

   TString pattern(fCurrentClass->GetName());
   // take unscoped version
   Ssiz_t posLastScope = kNPOS;
   while ((posLastScope = pattern.Index("::")) != kNPOS)
      pattern.Remove(0, posLastScope + 1);
   pattern += "::";

   TString declFileName;
   if (fHtml->GetDeclFileName(fCurrentClass, kTRUE, declFileName)) {
      LocateMethods(out, declFileName, kTRUE /*source info*/, useDocxxStyle,
                    kFALSE /*allowPureVirtual*/, pattern, 0);
      Ssiz_t posGt = pattern.Index('>');
      if (posGt != kNPOS) {
         // template! Re-run with pattern '...<.*>::'
         Ssiz_t posLt = pattern.Index('<');
         if (posLt != kNPOS && posLt < posGt) {
            pattern.Replace(posLt + 1, posGt - posLt - 1, ".*");
            LocateMethods(out, declFileName, kTRUE /*source info*/, useDocxxStyle,
                    kFALSE /*allowPureVirtual*/, pattern, 0);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Given fCurrentClass, look for methods in its header file's
/// class declaration block, and extract documentation to out,
/// while beautifying the header file in parallel.

void TDocParser::LocateMethodsInHeaderClassDecl(std::ostream& out)
{
   TString declFileName;
   if (fHtml->GetDeclFileName(fCurrentClass, kTRUE, declFileName))
      LocateMethods(out, declFileName, kTRUE/*source info*/, kTRUE /*useDocxxStyle*/,
                    kTRUE /*allowPureVirtual*/, 0, ".h.html");
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the current line as a comment, handling directives and re-formatting
/// the comment: remove "/*", "*/", "//", similar characters surrounding lines,
/// etc.
///
/// Return kFALSE if the line is not a comment.

Bool_t TDocParser::ProcessComment()
{
   if (!fCommentAtBOL
      && !(fLineStripped[0] == '/'
         && (fLineStripped[1] == '/' || fLineStripped[1] == '*'))
      && !InContext(kComment) && !InContext(kDirective)) {
      fLineComment = "";
      return kFALSE;
   }

   //if (InContext(kDirective))
   //   return kTRUE; - NO! we might have a comment from a previous directive!

   // don't write out empty lines if the current directive is eating the line
   if (InContext(kDirective) && !fLineComment.Length())
      return kTRUE;

   TString commentLine(fLineComment.Strip());

   // remove all <span class="comment"> tags
   Bool_t mustDealWithCommentAtBOL = fCommentAtBOL; // whether we had a closing "*/"
   Ssiz_t posComment = kNPOS;
   if (!fCommentAtBOL)
      posComment = commentLine.Index("<span class=\"comment\">", 0, TString::kIgnoreCase);
   Ssiz_t posSpanEnd = commentLine.Index("</span>", posComment == kNPOS?0:posComment, TString::kIgnoreCase);
   while ((mustDealWithCommentAtBOL && posSpanEnd != kNPOS) || posComment != kNPOS) {
      Int_t spanLevel = 1;
      Ssiz_t posSpan = commentLine.Index("<span", posComment + 1, TString::kIgnoreCase);
      while (spanLevel > 1 || (posSpan != kNPOS && posSpan < posSpanEnd)) {
         // another span was opened, take the next </span>
         if (posSpan != kNPOS && posSpan < posSpanEnd) {
            ++spanLevel;
            posSpan = commentLine.Index("<span", posSpan + 1, TString::kIgnoreCase);
            // posSpanEnd doesn't change
            continue;
         } // else
         --spanLevel;
         // posSpan doesn't change
         posSpanEnd = commentLine.Index("</span>", posSpanEnd + 1, TString::kIgnoreCase);
      }
      if (posSpanEnd != kNPOS) {
          // only remove span if </span> if it exists (or we end up with unbalanced spans)
         commentLine.Remove(posSpanEnd, 7);
         if (posComment != kNPOS)
            commentLine.Remove(posComment, 22);
         else {
            mustDealWithCommentAtBOL = kFALSE;
            // now remove C comments
            posComment = 0;
         }
         posComment = commentLine.Index("<span class=\"comment\">", posComment, TString::kIgnoreCase);
      } else break;
   }
   if (posComment != kNPOS)
      commentLine.Remove(posComment, 22);

   // don't strip in C comments, do strip if opening:
   if (!InContext(kComment) || (InContext(kComment) & kCXXComment)
       || (fLineStripped[0] == '/' && fLineStripped[1] == '*'))
      Strip(commentLine);

   // look for start tag of class description
   if ((fClassDocState == kClassDoc_LookingNothingFound
      || fClassDocState == kClassDoc_LookingHaveSomething)
      && !fComment.Length()
      && fDocContext == kIgnore && commentLine.Contains(fClassDescrTag)) {
      fDocContext = kDocClass;
   }

   char start_or_end = 0;
   // remove leading /*, //
   if (commentLine.Length()>1 && commentLine[0] == '/'
       && (commentLine[1] == '/' || commentLine[1] == '*')) {
      start_or_end = commentLine[1];
      commentLine.Remove(0, 2);
   }
   // remove trailing */
   if (start_or_end != '/' && commentLine.Length()>1
       && commentLine[commentLine.Length() - 2] == '*'
       && commentLine[commentLine.Length() - 1] == '/') {
      start_or_end = commentLine[commentLine.Length() - 2];
      commentLine.Remove(commentLine.Length()-2);
   }

   // remove repeating characters from the end of the line
   if (start_or_end && commentLine.Length() > 3) {
      TString lineAllOneChar(commentLine.Strip());

      Ssiz_t len = lineAllOneChar.Length();
      if (len > 2) {
         Char_t c = lineAllOneChar[len - 1];
         if (c == lineAllOneChar[len - 2] && c == lineAllOneChar[len - 3]) {
            TString lineAllOneCharStripped = lineAllOneChar.Strip(TString::kTrailing, c);
            Strip(lineAllOneCharStripped);
            if (!lineAllOneCharStripped.Length()) {
               commentLine.Remove(0);

               // also a class doc signature: line consists of ////
               if ((fClassDocState == kClassDoc_LookingNothingFound
                  || fClassDocState == kClassDoc_LookingHaveSomething)
                  && !fComment.Length()
                   && fDocContext == kIgnore && start_or_end=='/') {
                  fDocContext = kDocClass;
               }
            }
         }
      }
   }

   // remove leading and trailing chars from e.g. // some doc //
   if (commentLine.Length() > 0 && start_or_end == commentLine[commentLine.Length() - 1])
      // we already removed it as part of // or / *; also remove the trailing
      commentLine = commentLine.Strip(TString::kTrailing, start_or_end);

   if (commentLine.Length() > 2 && Context() != kDirective)
      while (commentLine.Length() > 2
             && !IsWord(commentLine[0])
             && commentLine[0] == commentLine[commentLine.Length() - 1])
         commentLine = commentLine.Strip(TString::kBoth, commentLine[0]);

   // remove leading '/' if we had // or '*' if we had / *
   while (start_or_end && commentLine[0] == start_or_end)
      commentLine.Remove(0, 1);

   fComment += commentLine + "\n";

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// remove the top-most comment context that matches cxxcomment,

void TDocParser::RemoveCommentContext(Bool_t cxxcomment)
{
   UInt_t lookFor = kComment;
   if (cxxcomment) lookFor |= kCXXComment;
   std::list<UInt_t>::iterator iComment = fParseContext.end();
   for (std::list<UInt_t>::iterator iContext = fParseContext.begin();
      iContext != fParseContext.end(); ++ iContext)
      if (*iContext == lookFor) iComment =iContext;
   if (iComment != fParseContext.end())
      fParseContext.erase(iComment);
}

////////////////////////////////////////////////////////////////////////////////
/// strips ' ', tabs, and newlines from both sides of str

Bool_t TDocParser::Strip(TString& str)
{
   Bool_t changed = str[0] == ' ' || str[0] == '\t' || str[0] == '\n';
   changed |= str.Length()
      && (str[str.Length() - 1] == ' ' || str[str.Length() - 1] == '\t'
         || str[str.Length() - 1] == '\n');
   if (!changed) return kFALSE;
   Ssiz_t i = 0;
   while (str[i] == ' ' || str[i] == '\t' || str[i] == '\n')
      ++i;
   str.Remove(0,i);
   i = str.Length() - 1;
   while (i >= 0 && (str[i] == ' ' || str[i] == '\t' || str[i] == '\n'))
      --i;
   str.Remove(i + 1, str.Length());
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Write the class description depending (among others) on fClassDocState.

void TDocParser::WriteClassDoc(std::ostream& out, Bool_t first /*= kTRUE*/)
{
   if (fClassDocState == kClassDoc_LookingHaveSomething || fClassDocState == kClassDoc_LookingNothingFound) {
      TString& classDoc = first || !fLastClassDoc.Length() ? fFirstClassDoc : fLastClassDoc;
      static_cast<TClassDocOutput*>(fDocOutput)->WriteClassDescription(out, classDoc);
      fClassDocState = kClassDoc_Written;
   }

}

namespace {
   static void RemoveUnneededSpaces(TString& s) {
      // Remove spaces except between identifier characters.
      // Assumes s is stripped (does not start nor end with space).
      for (Ssiz_t i = 1; i < s.Length() - 1; ++i) {
         if (s[i] == ' ') {
            char p = s[i - 1];
            char n = s[i + 1];
            if (((isalnum(p) || p == '_') && (isalnum(n) || n == '_'))
                || (p == '>' && n == '>')) {
               // "id id" or "> >": keep space
            } else {
               while (isspace(s[i])) {
                  s.Remove(i, 1);
               }
            }
         }
      }
   }

   static void ParseParameters(TString& strippedParams, TList& paramArr) {
      // Extract a list of strings (the parameters without initializers) from
      // the signature.
      int nest = 0;
      bool init = false;
      bool quoted = false;
      Ssiz_t len = strippedParams.Length();
      TString arg;
      for (Ssiz_t i = 0; i < len; ++i) {
         switch (strippedParams[i]) {
         case '<': // fallthrough
         case '(': // fallthrough
         case '[': ++nest; break;
         case '>': // fallthrough
         case ')': // fallthrough
         case ']': --nest; break;
         case '=': init = true; break;
         case '\'': ++i; if (strippedParams[i] == '\\') ++i; ++i; continue;
         case '\\': ++i; continue; break;
         case '"': quoted = !quoted; break;
         case ',': {
            if (!quoted && !nest) {
               TString strippedArg(arg.Strip(TString::kBoth));
               paramArr.AddLast(new TObjString(strippedArg));
               init = false;
               arg.Remove(0);
               continue;
            }
         }
         }
         if (!init) {
            arg += strippedParams[i];
         }
      }
      TString strippedLastArg(arg.Strip(TString::kBoth));
      if (strippedLastArg.Length()) {
         paramArr.AddLast(new TObjString(strippedLastArg));
      }
   }

   void MatchOverloadSignatures(TCollection* candidates, TList* paramArr)
   {
      // Check type identity of candidate signatures. For each argument, check whether it
      // reduces the list of candidates to > 0 elements.
      TList suppressed;
      TIter iCandidate(candidates);
      int nparams = paramArr->GetSize();
      for (int iparam = 0; iparam < nparams && candidates->GetSize() > 1; ++iparam) {
         TString& srcArg = ((TObjString*)paramArr->At(iparam))->String();
         TString noParName(srcArg);
         while (noParName.Length()
                && (isalnum(noParName[noParName.Length() - 1]) || noParName[noParName.Length() - 1] == '_'))
            noParName.Remove(noParName.Length() - 1);
         noParName = noParName.Strip(TString::kTrailing);

         if (noParName.Length()) {
            RemoveUnneededSpaces(noParName);
         }
         RemoveUnneededSpaces(srcArg);
         // comparison:
         // 0: strcmp
         // 1: source's parameter has last identifier (parameter name?) removed
         // 2: candidate type name contained in source parameter
         for (int comparison = 0; comparison < 5; ++comparison) {
            if (comparison == 1 && noParName == srcArg)
               // there is no parameter name to ignore
               continue;
            suppressed.Clear();
            iCandidate.Reset();
            TDocMethodWrapper* method = 0;
            while ((method = (TDocMethodWrapper*) iCandidate())) {
               TMethodArg* methArg = (TMethodArg*) method->GetMethod()->GetListOfMethodArgs()->At(iparam);
               TString sMethArg = methArg->GetFullTypeName();
               RemoveUnneededSpaces(sMethArg);
               bool matches = false;
               switch (comparison) {
               case 0: matches = (srcArg == sMethArg); break;
               case 1: matches = (noParName == sMethArg); break;
               case 2: matches = srcArg.Contains(sMethArg) || sMethArg.Contains(srcArg); break;
               }
               if (!matches) {
                  suppressed.Add(method);
               }
            }
            if (suppressed.GetSize()
                && suppressed.GetSize() < candidates->GetSize()) {
               candidates->RemoveAll(&suppressed);
               break;
            }
            if (!suppressed.GetSize()) {
               // we have a match, no point in trying a looser matching
               break;
            }
         }
      }
      if (candidates->GetSize() > 1) {
         // use TDocMethodWrapper::kDocumented bit
         suppressed.Clear();
         iCandidate.Reset();
         TDocMethodWrapper* method = 0;
         while ((method = (TDocMethodWrapper*) iCandidate())) {
            if (method->TestBit(TDocMethodWrapper::kDocumented)) {
               suppressed.AddLast(method);
            }
         }
         if (suppressed.GetSize()
             && suppressed.GetSize() < candidates->GetSize()) {
            candidates->RemoveAll(&suppressed);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a method, forwarding to TClassDocOutput

void TDocParser::WriteMethod(std::ostream& out, TString& ret,
                             TString& name, TString& params, Bool_t isconst,
                             const char* filename, TString& anchor,
                             TString& codeOneLiner)
{
   // if we haven't found the class description until now it's too late.
   if (fClassDocState < kClassDoc_Written)
      WriteClassDoc(out);

   TString strippedParams(params);
   if (strippedParams[0] == '(') {
      strippedParams.Remove(0, 1);
      strippedParams.Remove(strippedParams.Length() - 1);
      strippedParams = strippedParams.Strip(TString::kBoth);
   }

   TList paramArr;
   paramArr.SetOwner();
   ParseParameters(strippedParams, paramArr);
   int nparams = paramArr.GetSize();

   // Collect overload candidates
   TList candidates;
   for (int access = 0; access < 3; ++access) {
      const TList* methList = fMethods[access].GetListForObject(name);
      if (!methList) continue;

      TIter nextMethod(methList);
      TDocMethodWrapper* method = 0;
      while ((method = (TDocMethodWrapper *) nextMethod())) {
         if (name == method->GetName()
             && isconst == ((method->GetMethod()->Property() & kIsConstMethod) > 0)
             && method->GetMethod()->GetListOfMethodArgs()->GetSize() == nparams) {
            candidates.Add(method);
         }
      }
   }

   if (nparams && candidates.GetSize() > 1) {
      MatchOverloadSignatures(&candidates, &paramArr);
   }

   TDocMethodWrapper* guessedMethod = 0;
   if (candidates.GetSize() == 1) {
      guessedMethod = (TDocMethodWrapper*) candidates.First();
      guessedMethod->SetBit(TDocMethodWrapper::kDocumented);
   }

   static_cast<TClassDocOutput*>(fDocOutput)->WriteMethod(out, ret, name, params, filename, anchor,
                                                           fComment, codeOneLiner, guessedMethod);

   DecrementMethodCount(name);
   ret.Remove(0);
   name.Remove(0);
   params.Remove(0);
   anchor.Remove(0);
   fComment.Remove(0);

   fDocContext = kIgnore;
}

////////////////////////////////////////////////////////////////////////////////
/// Write fLineSource to out.
/// Adjust relative paths first.

void TDocParser::WriteSourceLine(std::ostream& out)
{
   fDocOutput->AdjustSourcePath(fLineSource);
   out << fLineSource << std::endl;

}
