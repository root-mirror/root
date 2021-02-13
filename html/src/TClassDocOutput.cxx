// @(#)root/html:$Id$
// Author: Axel Naumann 2007-01-09

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TClassDocOutput.h"

#include "TBaseClass.h"
#include "TClassEdit.h"
#include "TDataMember.h"
#include "TMethodArg.h"
#include "TDataType.h"
#include "TDocInfo.h"
#include "TDocParser.h"
#include "TError.h"
#include "THtml.h"
#include "TMethod.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TVirtualPad.h"
#include "TVirtualMutex.h"
#include "Riostream.h"
#include <sstream>

//______________________________________________________________________________
//
// Write the documentation for a class or namespace. The documentation is
// parsed by TDocParser and then passed to TClassDocOutput to generate
// the class doc header, the class description, members overview, and method
// documentation. All generic output functionality is in TDocOutput; it is
// re-used in this derived class.
//
// You usually do not use this class yourself; it is invoked indirectly by
// THtml. Customization of the output should happen via the interfaces defined
// by THtml.
//______________________________________________________________________________


ClassImp(TClassDocOutput);

////////////////////////////////////////////////////////////////////////////////
/// Create an object given the invoking THtml object, and the TClass
/// object that we will generate output for.

TClassDocOutput::TClassDocOutput(THtml& html, TClass* cl, TList* typedefs):
   TDocOutput(html), fHierarchyLines(0), fCurrentClass(cl),
   fCurrentClassesTypedefs(typedefs), fParser(0)
{
   fParser = new TDocParser(*this, fCurrentClass);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor, deletes fParser

TClassDocOutput::~TClassDocOutput()
{
   delete fParser;
}

////////////////////////////////////////////////////////////////////////////////
/// Create HTML files for a single class.
///

void TClassDocOutput::Class2Html(Bool_t force)
{
   gROOT->GetListOfGlobals(kTRUE);

   // create a filename
   TString filename(fCurrentClass->GetName());
   NameSpace2FileName(filename);

   gSystem->PrependPathName(fHtml->GetOutputDir(), filename);

   filename += ".html";

   if (!force && !IsModified(fCurrentClass, kSource)
       && !IsModified(fCurrentClass, kDoc)) {
      Printf(fHtml->GetCounterFormat(), "-no change-", fHtml->GetCounter(), filename.Data());
      return;
   }

   // open class file
   std::ofstream classFile(filename);

   if (!classFile.good()) {
      Error("Make", "Can't open file '%s' !", filename.Data());
      return;
   }

   Printf(fHtml->GetCounterFormat(), "", fHtml->GetCounter(), filename.Data());

   // write a HTML header for the classFile file
   WriteHtmlHeader(classFile, fCurrentClass->GetName(), "", fCurrentClass);
   WriteClassDocHeader(classFile);

   // copy .h file to the Html output directory
   TString declf;
   if (fHtml->GetDeclFileName(fCurrentClass, kTRUE, declf))
      CopyHtmlFile(declf);

   // process a '.cxx' file
   fParser->Parse(classFile);

   // write classFile footer
   WriteHtmlFooter(classFile, "",
      fParser->GetSourceInfo(TDocParser::kInfoLastUpdate),
      fParser->GetSourceInfo(TDocParser::kInfoAuthor),
      fParser->GetSourceInfo(TDocParser::kInfoCopyright));
}

////////////////////////////////////////////////////////////////////////////////
/// Write the list of functions

void TClassDocOutput::ListFunctions(std::ostream& classFile)
{
   // loop to get a pointers to method names

   classFile << std::endl << "<div id=\"functions\">" << std::endl;
   TString mangled(fCurrentClass->GetName());
   NameSpace2FileName(mangled);
   classFile << "<h2><a id=\"" << mangled
      << ":Function_Members\"></a>Function Members (Methods)</h2>" << std::endl;

   const char* tab4nbsp="&nbsp;&nbsp;&nbsp;&nbsp;";
   TString declFile;
   fHtml->GetDeclFileName(fCurrentClass, kFALSE, declFile);
   if (fCurrentClass->Property() & kIsAbstract)
      classFile << "&nbsp;<br /><b>"
                << tab4nbsp << "This is an abstract class, constructors will not be documented.<br />" << std::endl
                << tab4nbsp << "Look at the <a href=\""
                << gSystem->BaseName(declFile)
                << "\">header</a> to check for available constructors.</b><br />" << std::endl;

   Int_t minAccess = 0;
   if (fHtml->IsNamespace(fCurrentClass))
      minAccess = TDocParser::kPublic;
   for (Int_t access = TDocParser::kPublic; access >= minAccess; --access) {

      const TList* methods = fParser->GetMethods((TDocParser::EAccess)access);
      if (methods->GetEntries() == 0)
         continue;

      classFile << "<div class=\"access\" ";
      const char* accessID [] = {"priv", "prot", "publ"};
      const char* accesstxt[] = {"private", "protected", "public"};

      classFile << "id=\"func" << accessID[access] << "\"><b>"
         << accesstxt[access] << ":</b>" << std::endl
         << "<table class=\"func\" id=\"tabfunc" << accessID[access] << "\" cellspacing=\"0\">" << std::endl;

      TIter iMethWrap(methods);
      TDocMethodWrapper *methWrap = 0;
      while ((methWrap = (TDocMethodWrapper*) iMethWrap())) {
         const TMethod* method = methWrap->GetMethod();

         // it's a c'tor - Cint stores the class name as return type
         Bool_t isctor = (!strcmp(method->GetName(), method->GetReturnTypeName()));
         // it's a d'tor - Cint stores "void" as return type
         Bool_t isdtor = (!isctor && method->GetName()[0] == '~');

         classFile << "<tr class=\"func";
         if (method->GetClass() != fCurrentClass)
            classFile << "inh";
         classFile << "\"><td class=\"funcret\">";
         if (kIsVirtual & method->Property()) {
            if (!isdtor)
               classFile << "virtual ";
            else
               classFile << " virtual";
         }

         if (kIsStatic & method->Property())
            classFile << "static ";

         if (!isctor && !isdtor)
            fParser->DecorateKeywords(classFile, method->GetReturnTypeName());

         TString mangledM(method->GetClass()->GetName());
         NameSpace2FileName(mangledM);
         classFile << "</td><td class=\"funcname\"><a class=\"funcname\" href=\"";
         if (method->GetClass() != fCurrentClass) {
            TString htmlFile;
            fHtml->GetHtmlFileName(method->GetClass(), htmlFile);
            classFile << htmlFile;
         }
         classFile << "#" << mangledM;
         classFile << ":";
         mangledM = method->GetName();
         NameSpace2FileName(mangledM);
         Int_t overloadIdx = methWrap->GetOverloadIdx();
         if (overloadIdx) {
            mangledM += "@";
            mangledM += overloadIdx;
         }
         classFile << mangledM << "\">";
         if (method->GetClass() != fCurrentClass) {
            classFile << "<span class=\"baseclass\">";
            ReplaceSpecialChars(classFile, method->GetClass()->GetName());
            classFile << "::</span>";
         }
         ReplaceSpecialChars(classFile, method->GetName());
         classFile << "</a>";

         fParser->DecorateKeywords(classFile, const_cast<TMethod*>(method)->GetSignature());
         bool propSignal = false;
         bool propMenu   = false;
         bool propToggle = false;
         bool propGetter = false;
         if (method->GetTitle()) {
            propSignal = (strstr(method->GetTitle(), "*SIGNAL*"));
            propMenu   = (strstr(method->GetTitle(), "*MENU*"));
            propToggle = (strstr(method->GetTitle(), "*TOGGLE*"));
            propGetter = (strstr(method->GetTitle(), "*GETTER"));
            if (propSignal || propMenu || propToggle || propGetter) {
               classFile << "<span class=\"funcprop\">";
               if (propSignal) classFile << "<abbr title=\"emits a signal\">SIGNAL</abbr> ";
               if (propMenu) classFile << "<abbr title=\"has a popup menu entry\">MENU</abbr> ";
               if (propToggle) classFile << "<abbr title=\"toggles a state\">TOGGLE</abbr> ";
               if (propGetter) {
                  TString getter(method->GetTitle());
                  Ssiz_t posGetter = getter.Index("*GETTER=");
                  getter.Remove(0, posGetter + 8);
                  classFile << "<abbr title=\"use " + getter + "() as getter\">GETTER</abbr> ";
               }
               classFile << "</span>";
            }
         }
         classFile << "</td></tr>" << std::endl;
      }
      classFile << std::endl << "</table></div>" << std::endl;
   }

   classFile << "</div>" << std::endl; // class="functions"
}

////////////////////////////////////////////////////////////////////////////////
/// Write the list of data members and enums

void  TClassDocOutput::ListDataMembers(std::ostream& classFile)
{
   // make a loop on data members
   Bool_t haveDataMembers = (fParser->GetDataMembers(TDocParser::kPrivate)->GetEntries() ||
                             fParser->GetDataMembers(TDocParser::kProtected)->GetEntries() ||
                             fParser->GetDataMembers(TDocParser::kPublic)->GetEntries() ||
                             fParser->GetEnums(TDocParser::kPublic)->GetEntries() ||
                             fParser->GetEnums(TDocParser::kProtected)->GetEntries() ||
                             fParser->GetEnums(TDocParser::kPrivate)->GetEntries());

   if (!haveDataMembers) return;

   classFile << std::endl << "<div id=\"datamembers\">" << std::endl;
   TString mangled(fCurrentClass->GetName());
   NameSpace2FileName(mangled);
   classFile << "<h2><a name=\"" << mangled
      << ":Data_Members\"></a>Data Members</h2>" << std::endl;

   for (Int_t access = 5; access >= 0 && !fHtml->IsNamespace(fCurrentClass); --access) {
      const TList* datamembers = 0;
      if (access > 2) datamembers = fParser->GetEnums((TDocParser::EAccess) (access - 3));
      else datamembers = fParser->GetDataMembers((TDocParser::EAccess) access);
      if (datamembers->GetEntries() == 0)
         continue;

      classFile << "<div class=\"access\" ";
      const char* what = "data";
      if (access > 2) what = "enum";
      const char* accessID [] = {"priv", "prot", "publ"};
      const char* accesstxt[] = {"private", "protected", "public"};

      classFile << "id=\"" << what << accessID[access%3] << "\"><b>"
         << accesstxt[access%3] << ":</b>" << std::endl
         << "<table class=\"data\" id=\"tab" << what << accessID[access%3] << "\" cellspacing=\"0\">" << std::endl;

      TIter iDM(datamembers);
      TDataMember *member = 0;
      TString prevEnumName;
      Bool_t prevIsInh = kTRUE;

      while ((member = (TDataMember*) iDM())) {
         Bool_t haveNewEnum = access > 2 && prevEnumName != member->GetTypeName();
         if (haveNewEnum) {
            if (prevEnumName.Length()) {
               classFile << "<tr class=\"data";
               if (prevIsInh)
                  classFile << "inh";
               classFile << "\"><td class=\"datatype\">};</td><td></td><td></td></tr>" << std::endl;
            }
            prevEnumName = member->GetTypeName();
         }

         classFile << "<tr class=\"data";
         prevIsInh = (member->GetClass() != fCurrentClass);
         if (prevIsInh)
            classFile << "inh";
         classFile << "\"><td class=\"datatype\">";
         if (haveNewEnum) {
            TString enumName(member->GetTypeName());
            TString myScope(fCurrentClass->GetName());
            myScope += "::";
            enumName.ReplaceAll(myScope, "");
            if (enumName.EndsWith("::"))
               enumName += "<i>[unnamed]</i>";
            Ssiz_t startClassName = 0;
            if (!enumName.BeginsWith("enum "))
               classFile << "enum ";
            else
               startClassName = 5;

            Ssiz_t endClassName = enumName.Last(':'); // need template handling here!
            if (endClassName != kNPOS && endClassName > 0 && enumName[endClassName - 1] == ':') {
               // TClass* cl = fHtml->GetClass(TString(enumName(startClassName, endClassName - startClassName - 1)));
               TSubString substr(enumName(startClassName, endClassName - startClassName + 1));
               // if (cl)
                  // ReferenceEntity(substr, cl);
               enumName.Insert(substr.Start() + substr.Length(), "</span>");
               enumName.Insert(substr.Start(), "<span class=\"baseclass\">");
            }
            classFile << enumName << " { ";
         } else
            if (access < 3) {
               if (member->Property() & kIsStatic)
                  classFile << "static ";
               std::string shortTypeName(fHtml->ShortType(member->GetFullTypeName()));
               fParser->DecorateKeywords(classFile, shortTypeName.c_str());
            }

         TString mangledM(member->GetClass()->GetName());
         NameSpace2FileName(mangledM);
         classFile << "</td><td class=\"dataname\"><a ";
         if (member->GetClass() != fCurrentClass) {
            classFile << "href=\"";
            TString htmlFile;
            fHtml->GetHtmlFileName(member->GetClass(), htmlFile);
            classFile << htmlFile << "#";
         } else
            classFile << "name=\"";
         classFile << mangledM;
         classFile << ":";
         mangledM = member->GetName();
         NameSpace2FileName(mangledM);
         classFile << mangledM << "\">";
         if (member->GetClass() == fCurrentClass)
            classFile << "</a>";
         if (access < 3 && member->GetClass() != fCurrentClass) {
            classFile << "<span class=\"baseclass\">";
            ReplaceSpecialChars(classFile, member->GetClass()->GetName());
            classFile << "::</span>";
         }
         ReplaceSpecialChars(classFile, member->GetName());

         // Add the dimensions to "array" members
         for (Int_t indx = 0; indx < member->GetArrayDim(); ++indx)
            if (member->GetMaxIndex(indx) <= 0)
               break;
            else
               classFile << "[" << member->GetMaxIndex(indx) << "]";

         if (member->GetClass() != fCurrentClass)
            classFile << "</a>";
         classFile << "</td>";
         if (member->GetTitle() && member->GetTitle()[0]) {
            classFile << "<td class=\"datadesc\">";
            ReplaceSpecialChars(classFile, member->GetTitle());
         } else classFile << "<td>";
         classFile << "</td></tr>" << std::endl;
      } // for members

      if (prevEnumName.Length()) {
         classFile << "<tr class=\"data";
         if (prevIsInh)
            classFile << "inh";
         classFile << "\"><td class=\"datatype\">};</td><td></td><td></td></tr>" << std::endl;
      }
      classFile << std::endl << "</table></div>" << std::endl;
   } // for access

   classFile << "</div>" << std::endl; // datamembers
}

////////////////////////////////////////////////////////////////////////////////
/// This function builds the class charts for one class in GraphViz/Dot format,
/// i.e. the inheritance diagram, the include dependencies, and the library
/// dependency.
///
/// Input: out      - output file stream

Bool_t TClassDocOutput::ClassDotCharts(std::ostream& out)
{
   if (!fHtml->HaveDot())
      return kFALSE;

   TString title(fCurrentClass->GetName());
   NameSpace2FileName(title);

   TString dir("inh");
   gSystem->PrependPathName(fHtml->GetOutputDir(), dir);
   gSystem->MakeDirectory(dir);

   dir = "inhmem";
   gSystem->PrependPathName(fHtml->GetOutputDir(), dir);
   gSystem->MakeDirectory(dir);

   dir = "incl";
   gSystem->PrependPathName(fHtml->GetOutputDir(), dir);
   gSystem->MakeDirectory(dir);

   dir = "lib";
   gSystem->PrependPathName(fHtml->GetOutputDir(), dir);
   gSystem->MakeDirectory(dir);

   TString filenameInh(title);
   gSystem->PrependPathName("inh", filenameInh);
   gSystem->PrependPathName(fHtml->GetOutputDir(), filenameInh);
   filenameInh += "_Inh";
   if (!CreateDotClassChartInh(filenameInh + ".dot") || !RunDot(filenameInh, &out))
      return kFALSE;

   TString filenameInhMem(title);
   gSystem->PrependPathName("inhmem", filenameInhMem);
   gSystem->PrependPathName(fHtml->GetOutputDir(), filenameInhMem);
   filenameInhMem += "_InhMem";
   if (CreateDotClassChartInhMem(filenameInhMem + ".dot"))
      RunDot(filenameInhMem, &out);

   TString filenameIncl(title);
   gSystem->PrependPathName("incl", filenameIncl);
   gSystem->PrependPathName(fHtml->GetOutputDir(), filenameIncl);
   filenameIncl += "_Incl";
   if (CreateDotClassChartIncl(filenameIncl + ".dot"))
      RunDot(filenameIncl, &out);

   TString filenameLib(title);
   gSystem->PrependPathName("lib", filenameLib);
   gSystem->PrependPathName(fHtml->GetOutputDir(), filenameLib);
   filenameLib += "_Lib";
   if (CreateDotClassChartLib(filenameLib + ".dot"))
      RunDot(filenameLib, &out);

   out << "<div class=\"tabs\">" << std::endl
       << "<a id=\"img" << title << "_Inh\" class=\"tabsel\" href=\"inh/" << title << "_Inh.png\" onclick=\"javascript:return SetImg('Charts','inh/" << title << "_Inh.png');\">Inheritance</a>" << std::endl
       << "<a id=\"img" << title << "_InhMem\" class=\"tab\" href=\"inhmem/" << title << "_InhMem.png\" onclick=\"javascript:return SetImg('Charts','inhmem/" << title << "_InhMem.png');\">Inherited Members</a>" << std::endl
       << "<a id=\"img" << title << "_Incl\" class=\"tab\" href=\"incl/" << title << "_Incl.png\" onclick=\"javascript:return SetImg('Charts','incl/" << title << "_Incl.png');\">Includes</a>" << std::endl
       << "<a id=\"img" << title << "_Lib\" class=\"tab\" href=\"lib/" << title << "_Lib.png\" onclick=\"javascript:return SetImg('Charts','lib/" << title << "_Lib.png');\">Libraries</a><br/>" << std::endl
       << "</div><div class=\"classcharts\"><div class=\"classchartswidth\"></div>" << std::endl
       << "<img id=\"Charts\" alt=\"Class Charts\" class=\"classcharts\" usemap=\"#Map" << title << "_Inh\" src=\"inh/" << title << "_Inh.png\"/></div>" << std::endl;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// This function builds the class tree for one class in HTML
/// (inherited and succeeding classes, called recursively)
///
///
/// Input: out      - output file stream
///        classPtr - pointer to the class
///        dir      - direction to traverse tree: up, down or both
///

void TClassDocOutput::ClassHtmlTree(std::ostream& out, TClass * classPtr,
                          ETraverse dir, int depth)
{
   if (dir == kBoth) {
      out << "<!--INHERITANCE TREE-->" << std::endl;

      // draw class tree into nested tables recursively
      out << "<table><tr><td width=\"10%\"></td><td width=\"70%\">"
          << "<a href=\"ClassHierarchy.html\">Inheritance Chart</a>:</td></tr>";
      out << "<tr class=\"inhtree\"><td width=\"10%\"></td><td width=\"70%\">";

      out << "<table class=\"inhtree\"><tr><td>" << std::endl;
      out << "<table width=\"100%\" border=\"0\" ";
      out << "cellpadding =\"0\" cellspacing=\"2\"><tr>" << std::endl;
   } else {
      out << "<table><tr>";
   }

   ////////////////////////////////////////////////////////
   // Loop up to mother classes
   if (dir == kUp || dir == kBoth) {

      // make a loop on base classes
      TBaseClass *inheritFrom;
      TIter nextBase(classPtr->GetListOfBases());

      UInt_t bgcolor=255-depth*8;
      Bool_t first = kTRUE;
      while ((inheritFrom = (TBaseClass *) nextBase())) {

         if (first) {
            out << "<td><table><tr>" << std::endl;
            first = kFALSE;
         } else
            out << "</tr><tr>" << std::endl;
         out << "<td bgcolor=\""
            << Form("#%02x%02x%02x", bgcolor, bgcolor, bgcolor)
            << "\" align=\"right\">" << std::endl;
         // get a class
         TClass *classInh = fHtml->GetClass((const char *) inheritFrom->GetName());
         if (classInh)
            ClassHtmlTree(out, classInh, kUp, depth+1);
         else
            out << "<tt>"
                << (const char *) inheritFrom->GetName()
                << "</tt>";
         out << "</td>"<< std::endl;
      }
      if (!first) {
         out << "</tr></table></td>" << std::endl; // put it in additional row in table
         out << "<td>&larr;</td>";
      }
   }

   out << "<td>" << std::endl; // put it in additional row in table
   ////////////////////////////////////////////////////////
   // Output Class Name

   const char *className = classPtr->GetName();
   TString htmlFile;
   fHtml->GetHtmlFileName(classPtr, htmlFile);
   TString anchor(className);
   NameSpace2FileName(anchor);

   if (dir == kUp) {
      if (htmlFile) {
         out << "<center><tt><a name=\"" << anchor;
         out << "\" href=\"" << htmlFile << "\">";
         ReplaceSpecialChars(out, className);
         out << "</a></tt></center>" << std::endl;
      } else
         ReplaceSpecialChars(out, className);
   }

   if (dir == kBoth) {
      if (htmlFile.Length()) {
         out << "<center><big><b><tt><a name=\"" << anchor;
         out << "\" href=\"" << htmlFile << "\">";
         ReplaceSpecialChars(out, className);
         out << "</a></tt></b></big></center>" << std::endl;
      } else
         ReplaceSpecialChars(out, className);
   }

   out << "</td>" << std::endl; // put it in additional row in table

   ////////////////////////////////////////////////////////
   // Loop down to child classes

   if (dir == kDown || dir == kBoth) {

      // 1. make a list of class names
      // 2. use DescendHierarchy

      out << "<td><table><tr>" << std::endl;
      fHierarchyLines = 0;
      DescendHierarchy(out,classPtr,10);

      out << "</tr></table>";
      if (dir==kBoth && fHierarchyLines>=10)
         out << "</td><td align=\"left\">&nbsp;<a href=\"ClassHierarchy.html\">[more...]</a>";
      out<<"</td>" << std::endl;

      // free allocated memory
   }

   out << "</tr></table>" << std::endl;
   if (dir == kBoth)
      out << "</td></tr></table></td></tr></table>"<<std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// It makes a graphical class tree
///
///
/// Input: psCanvas - pointer to the current canvas
///        classPtr - pointer to the class
///

void TClassDocOutput::ClassTree(TVirtualPad * psCanvas, Bool_t force)
{
   if (!psCanvas || !fCurrentClass)
      return;

   TString filename(fCurrentClass->GetName());
   NameSpace2FileName(filename);

   gSystem->PrependPathName(fHtml->GetOutputDir(), filename);


   filename += "_Tree.pdf";

   if (IsModified(fCurrentClass, kTree) || force) {
      // TCanvas already prints pdf being saved
      // Printf(fHtml->GetCounterFormat(), "", "", filename);
      fCurrentClass->Draw("same");
      Int_t saveErrorIgnoreLevel = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kWarning;
      psCanvas->SaveAs(filename);
      gErrorIgnoreLevel = saveErrorIgnoreLevel;
   } else
      Printf(fHtml->GetCounterFormat(), "-no change-", "", filename.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Build the class tree for one class in GraphViz/Dot format
///
///
/// Input: filename - output dot file incl. path

Bool_t TClassDocOutput::CreateDotClassChartInh(const char* filename)
{
   std::ofstream outdot(filename);
   outdot << "strict digraph G {" << std::endl
      << "rankdir=RL;" << std::endl
      << "ranksep=2;" << std::endl
      << "nodesep=0;" << std::endl
      << "size=\"8,10\";" << std::endl
      << "ratio=auto;" << std::endl
      << "margin=0;" << std::endl
      << "node [shape=plaintext,fontsize=40,width=4,height=0.75];" << std::endl
      << "\"" << fCurrentClass->GetName() << "\" [shape=ellipse];" << std::endl;

   std::stringstream ssDep;
   std::list<TClass*> writeBasesFor;
   writeBasesFor.push_back(fCurrentClass);
   Bool_t haveBases = fCurrentClass->GetListOfBases() &&
      fCurrentClass->GetListOfBases()->GetSize();
   if (haveBases) {
      outdot << "{" << std::endl;
      while (!writeBasesFor.empty()) {
         TClass* cl = writeBasesFor.front();
         writeBasesFor.pop_front();
         if (cl != fCurrentClass) {
            outdot << "  \"" << cl->GetName() << "\"";
            const char* htmlFileName = fHtml->GetHtmlFileName(cl->GetName());
            if (htmlFileName)
               outdot << " [URL=\"" << htmlFileName << "\"]";
            outdot << ";" << std::endl;
         }
         if (cl->GetListOfBases() && cl->GetListOfBases()->GetSize()) {
            ssDep << "  \"" << cl->GetName() << "\" -> {";
            TIter iBase(cl->GetListOfBases());
            TBaseClass* base = 0;
            while ((base = (TBaseClass*)iBase())) {
               ssDep << " \"" << base->GetName() << "\";";
               writeBasesFor.push_back(base->GetClassPointer());
            }
            ssDep << "}" << std::endl;
         }
      }
      outdot << "}" << std::endl; // cluster
   }

   std::map<TClass*, Int_t> derivesFromMe;
   std::map<TClass*, unsigned int> entriesPerDerived;
   std::set<TClass*> wroteNode;
   wroteNode.insert(fCurrentClass);
   static const unsigned int maxClassesPerDerived = 20;
   fHtml->GetDerivedClasses(fCurrentClass, derivesFromMe);
   outdot << "{" << std::endl;
   for (Int_t level = 1; kTRUE; ++level) {
      Bool_t levelExists = kFALSE;
      for (std::map<TClass*, Int_t>::iterator iDerived = derivesFromMe.begin();
         iDerived != derivesFromMe.end(); ++iDerived) {
         if (iDerived->second != level) continue;
         levelExists = kTRUE;
         TIter iBaseOfDerived(iDerived->first->GetListOfBases());
         TBaseClass* baseDerived = 0;
         Bool_t writeNode = kFALSE;
         TClass* writeAndMoreFor = 0;
         while ((baseDerived = (TBaseClass*) iBaseOfDerived())) {
            TClass* clBaseDerived = baseDerived->GetClassPointer();
            if (clBaseDerived->InheritsFrom(fCurrentClass)
               && wroteNode.find(clBaseDerived) != wroteNode.end()) {
               unsigned int& count = entriesPerDerived[clBaseDerived];
               if (count < maxClassesPerDerived) {
                  writeNode = kTRUE;
                  ssDep << "\"" << iDerived->first->GetName() << "\" -> \""
                     << clBaseDerived->GetName() << "\";" << std::endl;
                  ++count;
               } else if (count == maxClassesPerDerived) {
                  writeAndMoreFor = clBaseDerived;
                  ssDep << "\"...andmore" << clBaseDerived->GetName() << "\"-> \""
                     << clBaseDerived->GetName() << "\";" << std::endl;
                  ++count;
               }
            }
         }

         if (writeNode) {
            wroteNode.insert(iDerived->first);
            outdot << "  \"" << iDerived->first->GetName() << "\"";
            const char* htmlFileName = fHtml->GetHtmlFileName(iDerived->first->GetName());
            if (htmlFileName)
               outdot << " [URL=\"" << htmlFileName << "\"]";
            outdot << ";" << std::endl;
         } else if (writeAndMoreFor) {
               outdot << "  \"...andmore" << writeAndMoreFor->GetName()
                      << "\" [label=\"...and more\",fontname=\"Times-Italic\",fillcolor=lightgrey,style=filled];" << std::endl;
         }
      }
      if (!levelExists) break;
   }
   outdot << "}" << std::endl; // cluster

   outdot << ssDep.str();

   outdot << "}" << std::endl; // digraph

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Build the class tree of inherited members for one class in GraphViz/Dot format
///
/// Input: filename - output dot file incl. path

Bool_t TClassDocOutput::CreateDotClassChartInhMem(const char* filename) {
   std::ofstream outdot(filename);
   outdot << "strict digraph G {" << std::endl
      << "ratio=auto;" << std::endl
      << "rankdir=RL;" << std::endl
      << "compound=true;" << std::endl
      << "constraint=false;" << std::endl
      << "ranksep=0.1;" << std::endl
      << "nodesep=0;" << std::endl
      << "margin=0;" << std::endl;
   outdot << "  node [style=filled,width=0.7,height=0.15,fixedsize=true,shape=plaintext,fontsize=10];" << std::endl;

   std::stringstream ssDep;
   const int numColumns = 3;

   std::list<TClass*> writeBasesFor;
   writeBasesFor.push_back(fCurrentClass);
   while (!writeBasesFor.empty()) {
      TClass* cl = writeBasesFor.front();
      writeBasesFor.pop_front();

      const char* htmlFileName = fHtml->GetHtmlFileName(cl->GetName());

      outdot << "subgraph \"cluster" << cl->GetName() << "\" {" << std::endl
             << "  color=lightgray;" << std::endl
             << "  label=\"" << cl->GetName() << "\";" << std::endl;
      if (cl != fCurrentClass && htmlFileName)
         outdot << "  URL=\"" << htmlFileName << "\"" << std::endl;

      //Bool_t haveMembers = (cl->GetListOfDataMembers() && cl->GetListOfDataMembers()->GetSize());
      Bool_t haveFuncs = cl->GetListOfMethods() && cl->GetListOfMethods()->GetSize();

      // DATA MEMBERS
      {
         // make sure each member name is listed only once
         // that's useless for data members, but symmetric to what we have for methods
         std::map<std::string, TDataMember*> dmMap;

         {
            TIter iDM(cl->GetListOfDataMembers());
            TDataMember* dm = 0;
            while ((dm = (TDataMember*) iDM()))
               dmMap[dm->GetName()] = dm;
         }

         outdot << "subgraph \"clusterData0" << cl->GetName() << "\" {" << std::endl
                << "  color=white;" << std::endl
                << "  label=\"\";" << std::endl
                << "  \"clusterNode0" << cl->GetName() << "\" [height=0,width=0,style=invis];" << std::endl;
         TString prevColumnNode;
         Int_t pos = dmMap.size();
         Int_t column = 0;
         Int_t newColumnEvery = (pos + numColumns - 1) / numColumns;
         for (std::map<std::string, TDataMember*>::iterator iDM = dmMap.begin();
              iDM != dmMap.end(); ++iDM, --pos) {
            TDataMember* dm = iDM->second;
            TString nodeName(cl->GetName());
            nodeName += "::";
            nodeName += dm->GetName();
            if (iDM == dmMap.begin())
               prevColumnNode = nodeName;

            outdot << "\"" << nodeName << "\" [label=\""
                   << dm->GetName() << "\"";
            if (dm->Property() & kIsPrivate)
               outdot << ",color=\"#FFCCCC\"";
            else if (dm->Property() & kIsProtected)
               outdot << ",color=\"#FFFF77\"";
            else
               outdot << ",color=\"#CCFFCC\"";
            outdot << "];" << std::endl;
            if (pos % newColumnEvery == 1) {
               ++column;
               outdot << "};" << std::endl // end dataR
                      << "subgraph \"clusterData" << column << cl->GetName() << "\" {" << std::endl
                      << "  color=white;" << std::endl
                      << "  label=\"\";" << std::endl
                      << "  \"clusterNode" << column << cl->GetName() << "\" [height=0,width=0,style=invis];" << std::endl;
            } else if (iDM != dmMap.begin() && pos % newColumnEvery == 0) {
               ssDep << "\"" << prevColumnNode
                     << "\" -> \"" << nodeName << "\""<< " [style=invis,weight=100];" << std::endl;
               prevColumnNode = nodeName;
            }
         }

         while (column < numColumns - 1) {
            ++column;
            outdot << "  \"clusterNode" << column << cl->GetName() << "\" [height=0,width=0,style=invis];" << std::endl;
         }

         outdot << "};" << std::endl; // subgraph dataL/R
      } // DATA MEMBERS

      // FUNCTION MEMBERS
      if (haveFuncs) {
         // make sure each member name is listed only once
         std::map<std::string, TMethod*> methMap;

         {
            TIter iMeth(cl->GetListOfMethods());
            TMethod* meth = 0;
            while ((meth = (TMethod*) iMeth()))
               methMap[meth->GetName()] = meth;
         }

         outdot << "subgraph \"clusterFunc0" << cl->GetName() << "\" {" << std::endl
                << "  color=white;" << std::endl
                << "  label=\"\";" << std::endl
                << "  \"clusterNode0" << cl->GetName() << "\" [height=0,width=0,style=invis];" << std::endl;

         TString prevColumnNodeFunc;
         Int_t pos = methMap.size();
         Int_t column = 0;
         Int_t newColumnEvery = (pos + numColumns - 1) / numColumns;
         for (std::map<std::string, TMethod*>::iterator iMeth = methMap.begin();
            iMeth != methMap.end(); ++iMeth, --pos) {
            TMethod* meth = iMeth->second;
            TString nodeName(cl->GetName());
            nodeName += "::";
            nodeName += meth->GetName();
            if (iMeth == methMap.begin())
               prevColumnNodeFunc = nodeName;

            outdot << "\"" << nodeName << "\" [label=\"" << meth->GetName() << "\"";
            if (cl != fCurrentClass &&
               fCurrentClass->GetMethodAny(meth->GetName()))
               outdot << ",color=\"#777777\"";
            else if (meth->Property() & kIsPrivate)
               outdot << ",color=\"#FFCCCC\"";
            else if (meth->Property() & kIsProtected)
               outdot << ",color=\"#FFFF77\"";
            else
               outdot << ",color=\"#CCFFCC\"";
            outdot << "];" << std::endl;
            if (pos % newColumnEvery == 1) {
               ++column;
               outdot << "};" << std::endl // end funcR
                      << "subgraph \"clusterFunc" << column << cl->GetName() << "\" {" << std::endl
                      << "  color=white;" << std::endl
                      << "  label=\"\";" << std::endl;
            } else if (iMeth != methMap.begin() && pos % newColumnEvery == 0) {
               ssDep << "\"" << prevColumnNodeFunc
                     << "\" -> \"" << nodeName << "\""<< " [style=invis,weight=100];" << std::endl;
               prevColumnNodeFunc = nodeName;
            }
         }
         outdot << "};" << std::endl; // subgraph funcL/R
      }

      outdot << "}" << std::endl; // cluster class

      for (Int_t pos = 0; pos < numColumns - 1; ++pos)
         ssDep << "\"clusterNode" << pos << cl->GetName() << "\" -> \"clusterNode" << pos + 1 << cl->GetName() << "\" [style=invis];" << std::endl;

      if (cl->GetListOfBases() && cl->GetListOfBases()->GetSize()) {
         TIter iBase(cl->GetListOfBases());
         TBaseClass* base = 0;
         while ((base = (TBaseClass*)iBase())) {
            ssDep << "  \"clusterNode" << numColumns - 1 << cl->GetName() << "\" -> "
                  << " \"clusterNode0" << base->GetName() << "\" [ltail=\"cluster" << cl->GetName()
                  << "\",lhead=\"cluster" << base->GetName() << "\"";
            if (base != cl->GetListOfBases()->First())
               ssDep << ",weight=0";
            ssDep << "];" << std::endl;
            writeBasesFor.push_back(base->GetClassPointer());
         }
      }
   }

   outdot << ssDep.str();

   outdot << "}" << std::endl; // digraph

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Build the include dependency graph for one class in
/// GraphViz/Dot format
///
/// Input: filename - output dot file incl. path

Bool_t TClassDocOutput::CreateDotClassChartIncl(const char* filename) {
   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   std::map<std::string, std::string> filesToParse;
   std::list<std::string> listFilesToParse;
   TString declFileName;
   TString implFileName;
   fHtml->GetImplFileName(fCurrentClass, kFALSE, implFileName);
   if (fHtml->GetDeclFileName(fCurrentClass, kFALSE, declFileName)) {
      TString real;
      if (fHtml->GetDeclFileName(fCurrentClass, kTRUE, real)) {
         filesToParse[declFileName.Data()] = real.Data();
         listFilesToParse.push_back(declFileName.Data());
      }
   }
   /* do it only for the header
   if (implFileName && strlen(implFileName)) {
      char* real = gSystem->Which(fHtml->GetInputPath(), implFileName, kReadPermission);
      if (real) {
         filesToParse[implFileName] = real;
         listFilesToParse.push_back(implFileName);
         delete real;
      }
   }
   */

   std::ofstream outdot(filename);
   outdot << "strict digraph G {" << std::endl
      << "ratio=compress;" << std::endl
      << "rankdir=TB;" << std::endl
      << "concentrate=true;" << std::endl
      << "ranksep=0;" << std::endl
      << "nodesep=0;" << std::endl
      << "size=\"8,10\";" << std::endl
      << "node [fontsize=20,shape=plaintext];" << std::endl;

   for (std::list<std::string>::iterator iFile = listFilesToParse.begin();
      iFile != listFilesToParse.end(); ++iFile) {
      std::ifstream in(filesToParse[*iFile].c_str());
      std::string line;
      while (in && !in.eof()) {
         std::getline(in, line);
         size_t pos = 0;
         while (line[pos] == ' ' || line[pos] == '\t') ++pos;
         if (line[pos] != '#') continue;
         ++pos;
         while (line[pos] == ' ' || line[pos] == '\t') ++pos;
         if (line.compare(pos, 8, "include ") != 0) continue;
         pos += 8;
         while (line[pos] == ' ' || line[pos] == '\t') ++pos;
         if (line[pos] != '"' && line[pos] != '<')
            continue;
         char delim = line[pos];
         if (delim == '<') delim = '>';
         ++pos;
         line.erase(0, pos);
         pos = 0;
         pos = line.find(delim);
         if (pos == std::string::npos) continue;
         line.erase(pos);
         if (filesToParse.find(line) == filesToParse.end()) {
            TString sysfilename;
            if (!GetHtml()->GetPathDefinition().GetFileNameFromInclude(line.c_str(), sysfilename))
               continue;
            listFilesToParse.push_back(line);
            filesToParse[line] = sysfilename;
            if (*iFile == implFileName.Data() || *iFile == declFileName.Data())
               outdot << "\"" << *iFile << "\" [style=filled,fillcolor=lightgray];" << std::endl;
         }
         outdot << "\"" << *iFile << "\" -> \"" << line << "\";" << std::endl;
      }
   }

   outdot << "}" << std::endl; // digraph

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Build the library dependency graph for one class in
/// GraphViz/Dot format
///
/// Input: filename - output dot file incl. path

Bool_t TClassDocOutput::CreateDotClassChartLib(const char* filename) {
   std::ofstream outdot(filename);
   outdot << "strict digraph G {" << std::endl
      << "ratio=auto;" << std::endl
      << "rankdir=RL;" << std::endl
      << "compound=true;" << std::endl
      << "constraint=false;" << std::endl
      << "ranksep=0.7;" << std::endl
      << "nodesep=0.3;" << std::endl
      << "size=\"8,8\";" << std::endl
      << "ratio=compress;" << std::endl;

   TString libs(fCurrentClass->GetSharedLibs());
   outdot << "\"All Libraries\" [URL=\"LibraryDependencies.html\",shape=box,rank=max,fillcolor=lightgray,style=filled];" << std::endl;

   if (libs.Length()) {
      TString firstLib(libs);
      Ssiz_t end = firstLib.Index(' ');
      if (end != kNPOS) {
         firstLib.Remove(end, firstLib.Length());
         libs.Remove(0, end + 1);
      } else libs = "";

      {
         Ssiz_t posExt = firstLib.First(".");
         if (posExt != kNPOS)
            firstLib.Remove(posExt, firstLib.Length());
      }

      outdot << "\"All Libraries\" -> \"" << firstLib << "\" [style=invis];" << std::endl;
      outdot << "\"" << firstLib << "\" -> {" << std::endl;

      if (firstLib != "libCore")
         libs += " libCore";
      if (firstLib != "libCint")
         libs += " libCint";
      TString thisLib;
      for (Ssiz_t pos = 0; pos < libs.Length(); ++pos)
         if (libs[pos] != ' ')
            thisLib += libs[pos];
         else if (thisLib.Length()) {
            Ssiz_t posExt = thisLib.First(".");
            if (posExt != kNPOS)
               thisLib.Remove(posExt, thisLib.Length());
            outdot << " \"" << thisLib << "\";";
            thisLib = "";
         }
      // remaining lib
      if (thisLib.Length()) {
         Ssiz_t posExt = thisLib.First(".");
         if (posExt != kNPOS)
            thisLib.Remove(posExt, thisLib.Length());
         outdot << " \"" << thisLib << "\";";
         thisLib = "";
      }
      outdot << "}" << std::endl; // dependencies
   } else
      outdot << "\"No rlibmap information available.\"" << std::endl;

   outdot << "}" << std::endl; // digraph

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the hierarchical class list part for the current class's
/// base classes. docFileName contains doc for fCurrentClass.
///

void TClassDocOutput::CreateClassHierarchy(std::ostream& out, const char* docFileName)
{
   // Find basic base classes
   TList *bases = fCurrentClass->GetListOfBases();
   if (!bases || bases->IsEmpty())
      return;

   out << "<hr />" << std::endl;

   out << "<table><tr><td><ul><li><tt>";
   if (docFileName) {
      out << "<a name=\"" << fCurrentClass->GetName() << "\" href=\""
          << docFileName << "\">";
      ReplaceSpecialChars(out, fCurrentClass->GetName());
      out << "</a>";
   } else {
      ReplaceSpecialChars(out, fCurrentClass->GetName());
   }

   // find derived classes
   out << "</tt></li></ul></td>";
   fHierarchyLines = 0;
   DescendHierarchy(out, fCurrentClass);

   out << "</tr></table>" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a hierarchical class list
/// The algorithm descends from the base classes and branches into
/// all derived classes. Mixing classes are displayed several times.
///
///

Bool_t TClassDocOutput::CreateHierarchyDot()
{
   const char* title = "ClassHierarchy";
   TString filename(title);
   gSystem->PrependPathName(fHtml->GetOutputDir(), filename);

   // open out file
   std::ofstream dotout(filename + ".dot");

   if (!dotout.good()) {
      Error("CreateHierarchy", "Can't open file '%s.dot' !",
            filename.Data());
      return kFALSE;
   }

   dotout << "digraph G {" << std::endl
          << "ratio=auto;" << std::endl
          << "rankdir=RL;" << std::endl;

   // loop on all classes
   TClassDocInfo* cdi = 0;
   TIter iClass(fHtml->GetListOfClasses());
   while ((cdi = (TClassDocInfo*)iClass())) {

      TDictionary *dict = cdi->GetClass();
      TClass *cl = dynamic_cast<TClass*>(dict);
      if (cl == 0) {
         if (!dict)
            Warning("THtml::CreateHierarchy", "skipping class %s\n", cdi->GetName());
         continue;
      }

      // Find immediate base classes
      TList *bases = cl->GetListOfBases();
      if (bases && !bases->IsEmpty()) {
         dotout << "\"" << cdi->GetName() << "\" -> { ";
         TIter iBase(bases);
         TBaseClass* base = 0;
         while ((base = (TBaseClass*) iBase())) {
            // write out current class
            if (base != bases->First())
               dotout << "; ";
            dotout << "\"" << base->GetName() << "\"";
         }
         dotout << "};" << std::endl;
      } else
         // write out current class - no bases
         dotout << "\"" << cdi->GetName() << "\";" << std::endl;

   }

   dotout << "}";
   dotout.close();

   std::ofstream out(filename + ".html");
   if (!out.good()) {
      Error("CreateHierarchy", "Can't open file '%s.html' !",
            filename.Data());
      return kFALSE;
   }

   Printf(fHtml->GetCounterFormat(), "", fHtml->GetCounter(), (filename + ".html").Data());
   // write out header
   WriteHtmlHeader(out, "Class Hierarchy");
   out << "<h1>Class Hierarchy</h1>" << std::endl;

   WriteSearch(out);

   RunDot(filename, &out);

   out << "<img usemap=\"#Map" << title << "\" src=\"" << title << ".png\"/>" << std::endl;
   // write out footer
   WriteHtmlFooter(out);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Open a Class.cxx.html file, where Class is defined by classPtr, and .cxx.html by extension
/// It's created in fHtml->GetOutputDir()/src. If successful, the HTML header is written to out.

void TClassDocOutput::CreateSourceOutputStream(std::ostream& out, const char* extension,
                                     TString& sourceHtmlFileName)
{
   TString sourceHtmlDir("src");
   gSystem->PrependPathName(fHtml->GetOutputDir(), sourceHtmlDir);
   // create directory if necessary
   {
      R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

      if (gSystem->AccessPathName(sourceHtmlDir))
         gSystem->MakeDirectory(sourceHtmlDir);
   }
   sourceHtmlFileName = fCurrentClass->GetName();
   NameSpace2FileName(sourceHtmlFileName);
   gSystem->PrependPathName(sourceHtmlDir, sourceHtmlFileName);
   sourceHtmlFileName += extension;
   dynamic_cast<std::ofstream&>(out).open(sourceHtmlFileName);
   if (!out) {
      Warning("LocateMethodsInSource", "Can't open beautified source file '%s' for writing!",
         sourceHtmlFileName.Data());
      sourceHtmlFileName.Remove(0);
      return;
   }

   // write a HTML header
   TString title(fCurrentClass->GetName());
   title += " - source file";
   WriteHtmlHeader(out, title, "../", fCurrentClass);
   out << "<div id=\"codeAndLineNumbers\"><pre class=\"listing\">" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Descend hierarchy recursively
/// loop over all classes and look for classes with base class basePtr

void TClassDocOutput::DescendHierarchy(std::ostream& out, TClass* basePtr, Int_t maxLines, Int_t depth)
{
   if (maxLines)
      if (fHierarchyLines >= maxLines) {
         out << "<td></td>" << std::endl;
         return;
      }

   UInt_t numClasses = 0;

   TClassDocInfo* cdi = 0;
   TIter iClass(fHtml->GetListOfClasses());
   while ((cdi = (TClassDocInfo*)iClass()) && (!maxLines || fHierarchyLines<maxLines)) {

      TClass *classPtr = dynamic_cast<TClass*>(cdi->GetClass());
      if (!classPtr) continue;

      // find base classes with same name as basePtr
      TList* bases=classPtr->GetListOfBases();
      if (!bases) continue;

      TBaseClass *inheritFrom=(TBaseClass*)bases->FindObject(basePtr->GetName());
      if (!inheritFrom) continue;

      if (!numClasses)
         out << "<td>&larr;</td><td><table><tr>" << std::endl;
      else
         out << "</tr><tr>"<<std::endl;
      fHierarchyLines++;
      numClasses++;
      UInt_t bgcolor=255-depth*8;
      out << "<td bgcolor=\""
          << Form("#%02x%02x%02x", bgcolor, bgcolor, bgcolor)
          << "\">";
      out << "<table><tr><td>" << std::endl;

      TString htmlFile(cdi->GetHtmlFileName());
      if (htmlFile.Length()) {
         out << "<center><tt><a name=\"" << cdi->GetName() << "\" href=\""
             << htmlFile << "\">";
         ReplaceSpecialChars(out, cdi->GetName());
         out << "</a></tt></center>";
      } else {
         ReplaceSpecialChars(out, cdi->GetName());
      }
      // write title
      // commented out for now because it reduces overview
      /*
        len = strlen(classNames[i]);
        for (Int_t w = 0; w < (maxLen - len + 2); w++)
        out << ".";
        out << " ";

        out << "<a name=\"Title:";
        out << classPtr->GetName();
        out << "\">";
        ReplaceSpecialChars(out, classPtr->GetTitle());
        out << "</a></tt>" << std::endl;
      */

      out << "</td>" << std::endl;
      DescendHierarchy(out,classPtr,maxLines, depth+1);
      out << "</tr></table></td>" << std::endl;

   }  // loop over all classes
   if (numClasses)
      out << "</tr></table></td>" << std::endl;
   else
      out << "<td></td>" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an output file with a graphical representation of the class
/// inheritance. If force, replace existing output file.
/// This routine does nothing if fHtml->HaveDot() is true - use
/// ClassDotCharts() instead!

void TClassDocOutput::MakeTree(Bool_t force /*= kFALSE*/)
{
   // class tree only if no dot, otherwise it's part of charts
   if (!fCurrentClass || fHtml->HaveDot())
      return;

   TString htmlFile;
   fHtml->GetHtmlFileName(fCurrentClass, htmlFile);
   if (htmlFile.Length()
       && (htmlFile.BeginsWith("http://")
           || htmlFile.BeginsWith("https://")
           || gSystem->IsAbsoluteFileName(htmlFile))
       ) {
      htmlFile.Remove(0);
   }

   if (!htmlFile.Length()) {
      TString what(fCurrentClass->GetName());
      what += " (source not found)";
      Printf(fHtml->GetCounterFormat(), "-skipped-", "", what.Data());
      return;
   }

   R__LOCKGUARD(GetHtml()->GetMakeClassMutex());

   // Create a canvas without linking against GUI libs
   Bool_t wasBatch = gROOT->IsBatch();
   if (!wasBatch)
      gROOT->SetBatch();
   TVirtualPad *psCanvas = (TVirtualPad*)gROOT->ProcessLineFast("new TCanvas(\"R__THtml\",\"psCanvas\",0,0,1000,1200);");
   if (!wasBatch)
      gROOT->SetBatch(kFALSE);

   if (!psCanvas) {
      Error("MakeTree", "Cannot create a TCanvas!");
      return;
   }

   // make a class tree
   ClassTree(psCanvas, force);

   psCanvas->Close();
   delete psCanvas;
}

////////////////////////////////////////////////////////////////////////////////
/// Called by TDocParser::LocateMethods(), this hook writes out the class description
/// found by TDocParser. It's even called if none is found, i.e. if the first method
/// has occurred before a class description is found, so missing class descriptions
/// can be handled.
/// For HTML, its creates the description block, the list of functions and data
/// members, and the inheritance tree or, if Graphviz's dot is found, the class charts.

void TClassDocOutput::WriteClassDescription(std::ostream& out, const TString& description)
{
   // Class Description Title
   out << "<div class=\"dropshadow\"><div class=\"withshadow\">";
   TString anchor(fCurrentClass->GetName());
   NameSpace2FileName(anchor);
   out << "<h1><a name=\"" << anchor;
   out << ":description\"></a>";

   if (fHtml->IsNamespace(fCurrentClass))
      out << "namespace ";
   else
      out << "class ";
   ReplaceSpecialChars(out, fCurrentClass->GetName());


   // make a loop on base classes
   Bool_t first = kTRUE;
   TBaseClass *inheritFrom;
   TIter nextBase(fCurrentClass->GetListOfBases());

   while ((inheritFrom = (TBaseClass *) nextBase())) {
      if (first) {
         out << ": ";
         first = kFALSE;
      } else
         out << ", ";
      Long_t property = inheritFrom->Property();
      if (property & kIsPrivate)
         out << "private ";
      else if (property & kIsProtected)
         out << "protected ";
      else
         out << "public ";

      // get a class
      TClass *classInh = fHtml->GetClass(inheritFrom->GetName());

      TString htmlFile;
      fHtml->GetHtmlFileName(classInh, htmlFile);

      if (htmlFile.Length()) {
         // make a link to the base class
         out << "<a href=\"" << htmlFile << "\">";
         ReplaceSpecialChars(out, inheritFrom->GetName());
         out << "</a>";
      } else
         ReplaceSpecialChars(out, inheritFrom->GetName());
   }
   out << "</h1>" << std::endl;

   out << "<div class=\"classdescr\">" << std::endl;

   if (description.Length())
      out << "<pre>" << description << "</pre>";

   // typedefs pointing to this class:
   if (fCurrentClassesTypedefs && !fCurrentClassesTypedefs->IsEmpty()) {
      out << "<h4>This class is also known as (typedefs to this class)</h4>";
      TIter iTD(fCurrentClassesTypedefs);
      bool firsttd = true;
      TDataType* dt = 0;
      while ((dt = (TDataType*) iTD())) {
         if (!firsttd)
            out << ", ";
         else firsttd = false;
         fParser->DecorateKeywords(out, dt->GetName());
      }
   }

   out << "</div>" << std::endl
       << "</div></div>" << std::endl;

   ListFunctions(out);
   ListDataMembers(out);

   // create dot class charts or an html inheritance tree
   out << "<h2><a id=\"" << anchor
      << ":Class_Charts\"></a>Class Charts</h2>" << std::endl;
   if (!fHtml->IsNamespace(fCurrentClass))
      if (!ClassDotCharts(out))
         ClassHtmlTree(out, fCurrentClass);

   // header for the following function docs:
   out << "<h2>Function documentation</h2>" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Write out the introduction of a class description (shortcuts and links)

void TClassDocOutput::WriteClassDocHeader(std::ostream& classFile)
{
   classFile << "<a name=\"TopOfPage\"></a>" << std::endl;


   // show box with lib, include
   // needs to go first to allow title on the left
   TString sTitle(fCurrentClass->GetName());
   ReplaceSpecialChars(sTitle);
   if (fHtml->IsNamespace(fCurrentClass))
      sTitle.Prepend("namespace ");
   else
      sTitle.Prepend("class ");

   TString sInclude;
   TString sLib;
   const char* lib=fCurrentClass->GetSharedLibs();
   GetHtml()->GetPathDefinition().GetIncludeAs(fCurrentClass, sInclude);
   if (lib) {
      char* libDup=StrDup(lib);
      char* libDupSpace=strchr(libDup,' ');
      if (libDupSpace) *libDupSpace=0;
      char* libDupEnd=libDup+strlen(libDup);
      while (libDupEnd!=libDup)
         if (*(--libDupEnd)=='.') {
            *libDupEnd=0;
            break;
         }
      sLib = libDup;
      delete[] libDup;
   }
   classFile << "<script type=\"text/javascript\">WriteFollowPageBox('"
             << sTitle << "','" << sLib << "','" << sInclude << "');</script>" << std::endl;

   TString modulename;
   fHtml->GetModuleNameForClass(modulename, fCurrentClass);
   TModuleDocInfo* module = (TModuleDocInfo*) fHtml->GetListOfModules()->FindObject(modulename);
   WriteTopLinks(classFile, module, fCurrentClass->GetName(), kFALSE);

   classFile << "<div class=\"descrhead\"><div class=\"descrheadcontent\">" << std::endl // descrhead line 3
      << "<span class=\"descrtitle\">Source:</span>" << std::endl;

   // make a link to the '.cxx' file
   TString classFileName(fCurrentClass->GetName());
   NameSpace2FileName(classFileName);

   TString headerFileName;
   fHtml->GetDeclFileName(fCurrentClass, kFALSE, headerFileName);
   TString sourceFileName;
   fHtml->GetImplFileName(fCurrentClass, kFALSE, sourceFileName);
   if (headerFileName.Length())
      classFile << "<a class=\"descrheadentry\" href=\"src/" << classFileName
                << ".h.html\">header file</a>" << std::endl;
   else
      classFile << "<a class=\"descrheadentry\"> </a>" << std::endl;

   if (sourceFileName.Length())
      classFile << "<a class=\"descrheadentry\" href=\"src/" << classFileName
                << ".cxx.html\">source file</a>" << std::endl;
   else
      classFile << "<a class=\"descrheadentry\"> </a>" << std::endl;

   if (!fHtml->IsNamespace(fCurrentClass) && !fHtml->HaveDot()) {
      // make a link to the inheritance tree (postscript)
      classFile << "<a class=\"descrheadentry\" href=\"" << classFileName << "_Tree.pdf\"";
      classFile << ">inheritance tree (.pdf)</a> ";
   }

   const TString& viewCVSLink = GetHtml()->GetViewCVS();
   Bool_t mustReplace = viewCVSLink.Contains("%f");
   if (viewCVSLink.Length()) {
      if (headerFileName.Length()) {
         TString link(viewCVSLink);
         TString sHeader(headerFileName);
         if (GetHtml()->GetProductName() && !strcmp(GetHtml()->GetProductName(), "ROOT")) {
            Ssiz_t posInclude = sHeader.Index("/include/");
            if (posInclude != kNPOS) {
               // Cut off ".../include", i.e. keep leading '/'
               sHeader.Remove(0, posInclude + 8);
            } else {
               // no /include/; maybe /inc?
               posInclude = sHeader.Index("/inc/");
               if (posInclude != kNPOS) {
                  sHeader = "/";
                  sHeader += sInclude;
               }
            }
            if (sourceFileName && strstr(sourceFileName, "src")) {
               TString src(sourceFileName);
               src.Remove(src.Index("src"), src.Length());
               src += "inc";
               sHeader.Prepend(src);
            } else {
               TString src(fCurrentClass->GetSharedLibs());
               Ssiz_t posEndLib = src.Index(' ');
               if (posEndLib != kNPOS)
                  src.Remove(posEndLib, src.Length());
               if (src.BeginsWith("lib"))
                  src.Remove(0, 3);
               posEndLib = src.Index('.');
               if (posEndLib != kNPOS)
                  src.Remove(posEndLib, src.Length());
               src.ToLower();
               src += "/inc";
               sHeader.Prepend(src);
            }
            if (sHeader.BeginsWith("tmva/inc/TMVA"))
               sHeader.Remove(8, 5);
         }
         if (mustReplace) link.ReplaceAll("%f", sHeader);
         else link += sHeader;
         classFile << "<a class=\"descrheadentry\" href=\"" << link << "\">viewVC header</a> ";
      } else
         classFile << "<a class=\"descrheadentry\"> </a> ";
      if (sourceFileName.Length()) {
         TString link(viewCVSLink);
         if (mustReplace) link.ReplaceAll("%f", sourceFileName);
         else link += sourceFileName;
         classFile << "<a class=\"descrheadentry\" href=\"" << link << "\">viewVC source</a> ";
      } else
         classFile << "<a class=\"descrheadentry\"> </a> ";
   }

   TString currClassNameMangled(fCurrentClass->GetName());
   NameSpace2FileName(currClassNameMangled);

   TString wikiLink = GetHtml()->GetWikiURL();
   if (wikiLink.Length()) {
      if (wikiLink.Contains("%c")) wikiLink.ReplaceAll("%c", currClassNameMangled);
      else wikiLink += currClassNameMangled;
      classFile << "<a class=\"descrheadentry\" href=\"" << wikiLink << "\">wiki</a> ";
   }

   classFile << std::endl << "</div></div>" << std::endl; // descrhead line 3

   classFile << "<div class=\"descrhead\"><div class=\"descrheadcontent\">" << std::endl // descrhead line 4
      << "<span class=\"descrtitle\">Sections:</span>" << std::endl
      << "<a class=\"descrheadentry\" href=\"#" << currClassNameMangled;
   if (fHtml->IsNamespace(fCurrentClass))
      classFile << ":description\">namespace description</a> ";
   else
      classFile << ":description\">class description</a> ";
   classFile << std::endl
      << "<a class=\"descrheadentry\" href=\"#" << currClassNameMangled << ":Function_Members\">function members</a>" << std::endl
      << "<a class=\"descrheadentry\" href=\"#" << currClassNameMangled << ":Data_Members\">data members</a>" << std::endl
      << "<a class=\"descrheadentry\" href=\"#" << currClassNameMangled << ":Class_Charts\">class charts</a>" << std::endl
      << "</div></div>" << std::endl // descrhead line 4
      << "</div>" << std::endl; // toplinks, from TDocOutput::WriteTopLinks

   WriteLocation(classFile, module, fCurrentClass->GetName());
}


////////////////////////////////////////////////////////////////////////////////
/// Write method name with return type ret and parameters param to out.
/// Build a link using file and anchor. Cooment it with comment, and
/// show the code codeOneLiner (set if the func consists of only one line
/// of code, immediately surrounded by "{","}"). Also updates fMethodNames's
/// count of method names.

void TClassDocOutput::WriteMethod(std::ostream& out, TString& ret,
                                  TString& name, TString& params,
                                  const char* filename, TString& anchor,
                                  TString& comment, TString& codeOneLiner,
                                  TDocMethodWrapper* guessedMethod)
{
   fParser->DecorateKeywords(ret);
   out << "<div class=\"funcdoc\"><span class=\"funcname\">"
      << ret << " <a class=\"funcname\" name=\"";
   TString mangled(fCurrentClass->GetName());
   NameSpace2FileName(mangled);
   out << mangled << ":";
   mangled = name;
   NameSpace2FileName(mangled);
   if (guessedMethod && guessedMethod->GetOverloadIdx()) {
      mangled += "@";
      mangled += guessedMethod->GetOverloadIdx();
   }
   out << mangled << "\" href=\"src/" << filename;
   if (anchor.Length())
      out << "#" << anchor;
   out << "\">";
   ReplaceSpecialChars(out, name);
   out << "</a>";
   if (guessedMethod) {
      out << "(";
      TMethodArg* arg;
      TIter iParam(guessedMethod->GetMethod()->GetListOfMethodArgs());
      Bool_t first = kTRUE;
      while ((arg = (TMethodArg*) iParam())) {
         if (!first) out << ", ";
         else first = kFALSE;
         TString paramGuessed(arg->GetFullTypeName());
         paramGuessed += " ";
         paramGuessed += arg->GetName();
         if (arg->GetDefault() && strlen(arg->GetDefault())) {
            paramGuessed += " = ";
            paramGuessed += arg->GetDefault();
         }
         fParser->DecorateKeywords(paramGuessed);
         out << paramGuessed;
      }
      out << ")";
      if (guessedMethod->GetMethod()->Property() & kIsConstMethod)
         out << " const";
   } else {
      fParser->DecorateKeywords(params);
      out << params;
   }
   out << "</span><br />" << std::endl;

   if (comment.Length())
      out << "<div class=\"funccomm\"><pre>" << comment << "</pre></div>" << std::endl;

   if (codeOneLiner.Length()) {
      out << std::endl << "<div class=\"code\"><code class=\"inlinecode\">"
          << codeOneLiner << "</code></div>" << std::endl
          << "<div style=\"clear:both;\"></div>" << std::endl;
      codeOneLiner.Remove(0);
   }
   out << "</div>" << std::endl;
}



