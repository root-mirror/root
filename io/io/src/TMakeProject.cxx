// @(#)root/io:$Id$
// Author: Rene Brun   12/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*
\class TMakeProject TMakeProject.cxx
\ingroup IO

Helper class implementing the TFile::MakeProject.
**/

#include <ctype.h>
#include "TMakeProject.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TList.h"
#include "TROOT.h"
#include "TMD5.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TError.h"

////////////////////////////////////////////////////////////////////////////////
/// Add an include statement, if it has not already been added.

void TMakeProject::AddUniqueStatement(FILE *fp, const char *statement, char *inclist)
{
   if (!strstr(inclist, statement)) {
      if (strlen(inclist)+strlen(statement) >= 50000) {
         Fatal("AddUniqueStatement","inclist too short need %u instead of 500000",UInt_t(strlen(inclist)+strlen(statement)));
      }
      strcat(inclist, statement);
      fprintf(fp, "%s", statement);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add an include statement, if it has not already been added.

void TMakeProject::AddInclude(FILE *fp, const char *header, Bool_t system, char *inclist)
{
   TString what;
   if (system) {
      what.Form("#include <%s>\n", header);
   } else {
      what.Form("#include \"%s\"\n", header);
   }
   AddUniqueStatement(fp, what.Data(), inclist);
}

////////////////////////////////////////////////////////////////////////////////
/// Chop the name by replacing the ending (before a potential extension) with
/// a md5 summary of the name.

void TMakeProject::ChopFileName(TString &name, Int_t limit)
{
   Ssiz_t len = name.Length();
   Bool_t has_extension = name.EndsWith(".h");
   if (has_extension)
      len -= 2;
   if (len >= limit) {
      if (has_extension) {
         name.Remove(name.Length()-2);
      }
      TMD5 md;
      md.Update((const UChar_t*)name.Data(),name.Length());
      md.Final();
      name.Remove( limit - 32 - 5); // Chop the part longer than 255 and keep space for the md5 and leave space for an extension
      name.Append( md.AsString() );
      if (has_extension) {
         name.Append( ".h" );
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Return the header name containing the description of name.

TString TMakeProject::GetHeaderName(const char *in_name, const TList *extrainfos, Bool_t includeNested)
{
   TString result;
   std::string strname( TClassEdit::GetLong64_Name( in_name ) );
   const char *name = strname.c_str();
   Int_t len = strlen(name);
   Int_t nest = 0;
   for (Int_t i = 0; i < len; ++i) {
      switch (name[i]) {
         case '<':
            ++nest;
            result.Append('_');
            break;
         case '>':
            --nest;
            result.Append('_');
            break;
         case ':':
            if (nest == 0 && name[i+1] == ':') {
               TString nsname(name, i);
               TClass *cl = gROOT->GetClass(nsname);
               Bool_t definedInParent = !includeNested && cl && (cl->Size() != 0 || (cl->Size()==0 && !cl->HasInterpreterInfo() /*empty 'base' class on file*/));
               if (!definedInParent && cl==0 && extrainfos!=0) {
                  TStreamerInfo *clinfo = (TStreamerInfo*)extrainfos->FindObject(nsname);
                  if (clinfo && clinfo->GetClassVersion() == -5) {
                     definedInParent = kTRUE;
                  }
               }
               if (definedInParent) {
                  // The requested class is actually nested inside
                  // the class whose name we already 'copied' to
                  // result.  The declaration will be in the same
                  // header file as the outer class.
                  if (strcmp(name + strlen(name) - 2, ".h") == 0) {
                     result.Append(".h");
                  }
                  ChopFileName(result,127);
                  return result;
               }
#ifndef WIN32
            }
            result.Append('_');
#else
            }
            if (name[i+1] == '/') {
               // don't replace the colon ':' in the case it's part of the drive name in a Windows path
               // e.g. don't convert from "c:/root/..." to "c_/root/..."
               result.Append(name[i]);
            } else {
               result.Append('_');
            }
#endif
            break;
         case ',':
         case '*':
         case '[':
         case ']':
         case ' ':
         case '(':
         case ')':
            result.Append('_');
            break;
         case '/':
         case '\\':
         default:
            result.Append(name[i]);
      }
   }
   ChopFileName(result,127);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Write the start of the class (forward) declaration.
/// If 'implementEmptyClass' is 3 then never add a #pragma

UInt_t TMakeProject::GenerateClassPrefix(FILE *fp, const char *clname, Bool_t top, TString &protoname,
      UInt_t *numberOfClasses, Int_t implementEmptyClass, Bool_t needGenericTemplate)
{
   // First open the namespace (if any)
   Int_t numberOfNamespaces = 0;
   const char *fullname = clname;

   Bool_t istemplate = kFALSE;
   if (strchr(clname, ':')) {
      // We might have a namespace in front of the classname.
      Int_t len = strlen(clname);
      const char *name = clname;
      UInt_t nest = 0;
      for (Int_t cur = 0; cur < len; ++cur) {
         switch (clname[cur]) {
            case '<':
               ++nest;
               istemplate = kTRUE;
               break;
            case '>':
               if (nest) --nest;
               break;
            case ':': {
                  if (nest == 0 && clname[cur+1] == ':') {
                     // We have a scope
                     TString nsname(clname, cur);
                     TClass *cl = gROOT->GetClass(nsname);
                     if (top) {
                        if (cl == 0 || (cl && cl->Size() == 0)) {
                           TString last(name, cur - (name - clname));
                           if ((numberOfClasses == 0 || *numberOfClasses == 0) && strchr(last.Data(), '<') == 0) {
                              fprintf(fp, "namespace %s {\n", last.Data());
                              ++numberOfNamespaces;
                           } else {
                              TString headername(GetHeaderName(last,0));
                              fprintf(fp, "#ifndef %s_h\n", headername.Data());
                              fprintf(fp, "#define %s_h\n", headername.Data());
                              GenerateClassPrefix(fp, last.Data(), top, protoname, 0);
                              fprintf(fp, "{\n");
                              fprintf(fp, "public:\n");
                              if (numberOfClasses) ++(*numberOfClasses);
                              istemplate = kFALSE;
                           }
                           name = clname + cur + 2;
                        }
                     } else {
                        istemplate = kFALSE;
                        name = clname + cur + 2;
                     }
                  }
                  break;
               }
         }
      }
      clname = name;
   } else {
      istemplate = strstr(clname, "<") != 0;
   }

   protoname = clname;

   if (implementEmptyClass==1) {
      TString headername(GetHeaderName(fullname,0));
      fprintf(fp, "#ifndef %s_h\n", headername.Data());
      fprintf(fp, "#define %s_h\n", headername.Data());
   }
   if (istemplate) {
      std::vector<const char*> argtype;

      Ssiz_t pos = protoname.First('<');
      UInt_t nparam = 1;
      if (pos != kNPOS) {
         if (isdigit(protoname[pos+1])) {
            argtype.push_back("int");
         } else {
            argtype.push_back("typename");
         }
         UInt_t nest = 0;
         for (Ssiz_t i = pos; i < protoname.Length(); ++i) {
            switch (protoname[i]) {
               case '<':
                  ++nest;
                  break;
               case '>':
                  if (nest) --nest;
                  break;
               case ',':
                  if (nest == 1) {
                     if (isdigit(protoname[i+1])) {
                        argtype.push_back("int");
                     } else {
                        argtype.push_back("typename");
                     }
                     ++nparam;
                  }
                  break;
            }
         }
         protoname.Remove(pos);
      }

      // Forward declaration of template.
      fprintf(fp, "template <");
      for (UInt_t p = 0; p < nparam; ++p) {
         if (p >= argtype.size() ) {
            fprintf(fp, "/* missing */ T%d", p);
         } else {
            fprintf(fp, "%s T%d", argtype[p], p);
         }
         if (p != (nparam - 1)) fprintf(fp, ", ");
      }
      if (needGenericTemplate) {
         fprintf(fp, "> class %s", protoname.Data());
      } else {
         fprintf(fp, "> class %s;\n", protoname.Data());
         fprintf(fp, "template <> ");
      }
   }

   if (implementEmptyClass) {
      if (istemplate) {
         if (!needGenericTemplate) {
            fprintf(fp, "class %s", clname);
         }
         fprintf(fp, " {\n");
         if (numberOfClasses) ++(*numberOfClasses);
         fprintf(fp, "public:\n");
         fprintf(fp, "operator int() { return 0; };\n");
      } else {
         fprintf(fp, "enum %s { kDefault_%s };\n", clname, clname);
         // The nesting space of this class may not be #pragma declared (and without it
         // the dictionary is broken), so for now skip those
         if (implementEmptyClass==1) {
            if (strchr(fullname, ':') == 0) {
               // yes this is too aggressive, this needs to be fixed properly by moving the #pragma out of band.
               fprintf(fp, "#ifdef __MAKECINT__\n#pragma link C++ class %s+;\n#endif\n", fullname);
            }
            fprintf(fp, "#endif\n");
         }
      }
   } else {
      if (!(istemplate && needGenericTemplate)) {
          fprintf(fp, "class %s", clname);
      }
   }
   return numberOfNamespaces;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an empty StreamerInfo for the given type (no recursion) if it is not
/// not known in the list of class.
///
/// If the type itself is a template,
/// we mark it with version 1 (a class) otherwise we mark it as version -3 (an enum).

void TMakeProject::GenerateMissingStreamerInfo(TList *extrainfos, const char *clname, Bool_t iscope)
{
   if (!TClassEdit::IsStdClass(clname) && !TClass::GetClass(clname) && gROOT->GetType(clname) == 0) {

      TStreamerInfo *info = (TStreamerInfo*)extrainfos->FindObject(clname);
      if (!info) {
         // The class does not exist, let's create it
         TStreamerInfo *newinfo = new TStreamerInfo();
         newinfo->SetName(clname);
         if (clname[strlen(clname)-1]=='>') {
            newinfo->SetTitle("Generated by MakeProject as an empty class template instantiation");
            newinfo->SetClassVersion(1);
         } else if (iscope) {
            newinfo->SetTitle("Generated by MakeProject as a namespace");
            newinfo->SetClassVersion(-4 /*namespace*/);
         } else {
            newinfo->SetTitle("Generated by MakeProject as an enum");
            newinfo->SetClassVersion(-3 /*enum*/);
         }
         extrainfos->Add(newinfo);
      } else {
         if (iscope) {
            if (info->GetClassVersion() == -3) {
               // This was marked as an enum but is also used as a scope,
               // so it was actually a class.
               info->SetTitle("Generated by MakeProject as an empty class");
               info->SetClassVersion(-5 /*class*/);
            }
         } else {
            if (info->GetClassVersion() == -4) {
               // This was marked as a 'namespace' but it is also used as a template parameter,
               // so it was actually a class.
               info->SetTitle("Generated by MakeProject as an empty class");
               info->SetClassVersion(-5 /*class*/);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an empty StreamerInfo for types that are used in templates parameters
/// but are not known in the list of class.
///
/// If the type itself is a template, we mark it with version 1 (a class)
/// otherwise we mark it as version -3 (an enum).

void TMakeProject::GenerateMissingStreamerInfos(TList *extrainfos, const char *clname)
{
   UInt_t len = strlen(clname);
   UInt_t nest = 0;
   UInt_t last = 0;
   //Bool_t istemplate = kFALSE; // mark whether the current right most entity is a class template.

   for (UInt_t i = 0; i < len; ++i) {
      switch (clname[i]) {
         case ':':
            if (nest == 0 && clname[i+1] == ':') {
               TString incName(clname, i);
               GenerateMissingStreamerInfo(extrainfos, incName.Data(), kTRUE);
               //istemplate = kFALSE;
            }
            break;
         case '<':
            ++nest;
            if (nest == 1) last = i + 1;
            break;
         case '>':
            if (nest == 0) return; // The name is not well formed, give up.
            --nest; /* intentional fall through to the next case */
         case ',':
            if ((clname[i] == ',' && nest == 1) || (clname[i] == '>' && nest == 0)) {
               TString incName(clname + last, i - last);
               incName = TClassEdit::ShortType(incName.Data(), TClassEdit::kDropTrailStar | TClassEdit::kLong64);
               if (clname[i] == '>' && nest == 1) incName.Append(">");

               if (isdigit(incName[0])) {
                  // Not a class name, nothing to do.
               } else {
                  GenerateMissingStreamerInfos(extrainfos,incName.Data());
               }
               last = i + 1;
            }
      }
   }
   GenerateMissingStreamerInfo(extrainfos,TClassEdit::ShortType(clname, TClassEdit::kDropTrailStar | TClassEdit::kLong64).c_str(),kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate an empty StreamerInfo for types that are used in templates parameters
/// but are not known in the list of class.
///
/// If the type itself is a template,
/// we mark it with version 1 (a class) otherwise we mark it as version -3 (an enum).

void TMakeProject::GenerateMissingStreamerInfos(TList *extrainfos, TStreamerElement *element)
{
   if (element->IsBase()) {
      TClass *elemCl = element->GetClassPointer();
      if (elemCl) GenerateMissingStreamerInfos(extrainfos,elemCl->GetName());
      else GenerateMissingStreamerInfos(extrainfos,element->GetTypeName());
   } else {
      GenerateMissingStreamerInfos(extrainfos,element->GetTypeName());
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Insert a (complete) forward declaration for the class 'clname'

UInt_t TMakeProject::GenerateForwardDeclaration(FILE *fp, const char *clname, char *inclist, Bool_t implementEmptyClass, Bool_t needGenericTemplate, const TList *extrainfos)
{
   UInt_t ninc = 0;

   if (strchr(clname, '<')) {
      ninc += GenerateIncludeForTemplate(fp, clname, inclist, kTRUE, extrainfos);
   }
   TString protoname;
   UInt_t numberOfClasses = 0;
   UInt_t numberOfNamespaces = GenerateClassPrefix(fp, clname, kTRUE, protoname, &numberOfClasses, implementEmptyClass, needGenericTemplate);

   if (!implementEmptyClass) fprintf(fp, ";\n");
   for (UInt_t i = 0;i < numberOfClasses;++i) {
      fprintf(fp, "}; // end of class.\n");
      fprintf(fp, "#endif\n");
   }
   for (UInt_t i = 0;i < numberOfNamespaces;++i) {
      fprintf(fp, "} // end of namespace.\n");
   }

   return ninc;
}

////////////////////////////////////////////////////////////////////////////////
/// Add to the header file, the \#include needed for the argument of
/// this template.

UInt_t TMakeProject::GenerateIncludeForTemplate(FILE *fp, const char *clname, char *inclist, Bool_t forward, const TList *extrainfos)
{
   UInt_t ninc = 0;
   UInt_t len = strlen(clname);
   UInt_t nest = 0;
   UInt_t last = 0;


   for (UInt_t i = 0; i < len; ++i) {
      switch (clname[i]) {
         case '<':
            ++nest;
            if (nest == 1) last = i + 1;
            break;
         case '>':
            if (nest==0) return ninc; // the name is not well formed, give up.
            --nest; /* intentional fall through to the next case */
         case ',':
            if ((clname[i] == ',' && nest == 1) || (clname[i] == '>' && nest == 0)) {
               TString incName(clname + last, i - last);
               incName = TClassEdit::ShortType(incName.Data(), TClassEdit::kDropTrailStar | TClassEdit::kLong64);
               if (clname[i] == '>' && nest == 1) incName.Append(">");
               Int_t stlType;
               if (isdigit(incName[0])) {
                  // Not a class name, nothing to do.
               } else if ((stlType = TClassEdit::IsSTLCont(incName))) {
                  const char *what = "";
                  switch (TMath::Abs(stlType))  {
                     case ROOT::kSTLvector:
                        what = "vector";
                        break;
                     case ROOT::kSTLlist:
                        what = "list";
                        break;
                     case ROOT::kSTLforwardlist:
                        what = "forward_list";
                        break;
                     case ROOT::kSTLdeque:
                        what = "deque";
                        break;
                     case ROOT::kSTLmap:
                     case ROOT::kSTLmultimap:
                        what = "map";
                        break;
                     case ROOT::kSTLunorderedmap:
                     case ROOT::kSTLunorderedmultimap:
                        what = "unordered_map";
                        break;
                     case ROOT::kSTLset:
                     case ROOT::kSTLmultiset:
                        what = "set";
                        break;
                     case ROOT::kSTLunorderedset:
                     case ROOT::kSTLunorderedmultiset:
                        what = "unordered_set";
                        break;
                     case ROOT::kSTLbitset:
                        what = "bitset";
                        break;
                     default:
                        what = "undetermined_stl_container";
                        break;
                  }
                  AddInclude(fp, what, kTRUE, inclist);
                  fprintf(fp, "namespace std {} using namespace std;\n");
                  ninc += GenerateIncludeForTemplate(fp, incName, inclist, forward, extrainfos);
               } else if (TClassEdit::IsStdPair(incName)) {
                  AddInclude(fp, "utility", kTRUE, inclist);
                  ninc += GenerateIncludeForTemplate(fp, incName, inclist, forward, extrainfos);
               } else if (strncmp(incName.Data(), "auto_ptr<", strlen("auto_ptr<")) == 0) {
                  AddInclude(fp, "memory", kTRUE, inclist);
                  ninc += GenerateIncludeForTemplate(fp, incName, inclist, forward, extrainfos);
               } else if (TClassEdit::IsStdClass(incName)) {
                  // Do nothing.
               } else {
                  TClass *cl = gROOT->GetClass(incName);
                  if (!forward && cl) {
                     if (cl->HasInterpreterInfo()) {
                        // We have the real dictionary for this class.

                        const char *include = cl->GetDeclFileName();
                        if (include && include[0]) {

                           if (strncmp(include, "include/", 8) == 0) {
                              include += 8;
                           }
                           if (strncmp(include, "include\\", 9) == 0) {
                              include += 9;
                           }
                           TMakeProject::AddInclude(fp, include, kFALSE, inclist);
                        }
                        GenerateIncludeForTemplate(fp, incName, inclist, forward, extrainfos);
                     } else {
                        incName = GetHeaderName(incName,extrainfos);
                        incName.Append(".h");
                        AddInclude(fp, incName, kFALSE, inclist);
                     }
                  } else if (incName.Length() && incName[0] != ' ' && gROOT->GetType(incName) == 0) {
                     Bool_t emptyclass = !cl;
                     if (emptyclass && extrainfos) {
                        TStreamerInfo *info = (TStreamerInfo*)extrainfos->FindObject(incName);
                        if (info && info->GetClassVersion() == -5) {
                           emptyclass = kFALSE;
                        }
                     }
                     GenerateForwardDeclaration(fp, incName, inclist, emptyclass, kFALSE, extrainfos);
                  }
               }
               last = i + 1;
            }
      }
   }

   Int_t stlType = TClassEdit::IsSTLCont(clname);
   if (stlType) {
      std::vector<std::string> inside;
      int nestedLoc;
      TClassEdit::GetSplit( clname, inside, nestedLoc, TClassEdit::kLong64 );
      Int_t stlkind =  TClassEdit::STLKind(inside[0]);
      TClass *key = TClass::GetClass(inside[1].c_str());
      if (key) {
         TString what;
         switch (stlkind)  {
            case ROOT::kSTLmap:
            case ROOT::kSTLmultimap: {
                  what = "pair<";
                  what += UpdateAssociativeToVector( inside[1].c_str() );
                  what += ",";
                  what += UpdateAssociativeToVector( inside[2].c_str() );
                  what += " >";
                  what.ReplaceAll("std::","");
                  // Only ask for it if needed.
                  TClass *paircl = TClass::GetClass(what.Data());
                  if (paircl == 0 || !paircl->HasInterpreterInfo()) {
                     AddUniqueStatement(fp, TString::Format("#ifdef __MAKECINT__\n#pragma link C++ class %s+;\n#endif\n", what.Data()), inclist);
                  }
                  break;
               }
         }
      }
   }

   return ninc;
}


////////////////////////////////////////////////////////////////////////////////
/// Add to the header file anything that need to appear after the class
/// declaration (this includes some \#pragma link).

void TMakeProject::GeneratePostDeclaration(FILE *fp, const TVirtualStreamerInfo *info, char *inclist)
{
   TIter next(info->GetElements());
   TStreamerElement *element;
   while( (element = (TStreamerElement*)next()) ) {
      Int_t stlType = TClassEdit::IsSTLCont(element->GetTypeName());
      if (stlType) {
         std::vector<std::string> inside;
         int nestedLoc;
         TClassEdit::GetSplit( element->GetTypeName(), inside, nestedLoc, TClassEdit::kLong64 );
         Int_t stlkind =  TClassEdit::STLKind(inside[0]);
         TClass *key = TClass::GetClass(inside[1].c_str());
         TString what;
         if (TClassEdit::IsStdPair(inside[1])) {
            what = inside[1].c_str();
         } else if (key) {
            switch (stlkind)  {
               case ROOT::kSTLmap:
               case ROOT::kSTLmultimap:
               {
                  // Already done (see GenerateIncludeForTemplate
                  break;
               }
               default:
                  break;
            }
         }
         if (what.Length()) {
            // Only ask for it if needed.
            TClass *paircl = TClass::GetClass(what.Data());
            if (paircl == 0 || !paircl->HasInterpreterInfo()) {
               AddUniqueStatement(fp, TString::Format("#ifdef __MAKECINT__\n#pragma link C++ class %s+;\n#endif\n",what.Data()), inclist);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If we have a map, multimap, set or multiset, plus unordered partners,
/// and the key is a class, we need to replace the
/// container by a vector since we don't have the
/// comparator function.
/// The 'name' is modified to return the change in the name,
/// if any.

static constexpr int str_length(const char* str)
{
    return *str ? 1 + str_length(str + 1) : 0;
}

TString TMakeProject::UpdateAssociativeToVector(const char *name)
{
   TString newname( name );

   constexpr auto auto_ptr_len = str_length("auto_ptr<");
   if (strncmp(name, "auto_ptr<", auto_ptr_len) == 0) {
      newname = "unique_ptr<";
      newname += (name + auto_ptr_len);
   } else if (strchr(name,'<')!=0) {
      std::vector<std::string> inside;
      int nestedLoc;
      unsigned int narg = TClassEdit::GetSplit( name, inside, nestedLoc, TClassEdit::kLong64 );

      Int_t stlkind =  TMath::Abs(TClassEdit::STLKind(inside[0]));

      for(unsigned int i = 1; i<narg; ++i) {
         inside[i] = UpdateAssociativeToVector( inside[i].c_str() );
      }

      if (nestedLoc) narg = nestedLoc;

      // Treat the trailing stars the same as nested loc (i.e. ends up, properly, tacking them up back at the end of the name)
      if (!inside[narg-1].empty() && inside[narg-1][0] == '*')
         narg = narg - 1;


      // Remove default allocator if any.
      static const char* allocPrefix = "std::allocator<";
      static const unsigned int allocPrefixLen (strlen(allocPrefix));
      switch (stlkind) {
         case ROOT::kSTLvector:
         case ROOT::kSTLlist:
         case ROOT::kSTLforwardlist:
         case ROOT::kSTLdeque:
            if (narg>2 && strncmp(inside[2].c_str(),allocPrefix,allocPrefixLen)==0) {
               --narg;
            }
            break;
         case ROOT::kSTLset:
         case ROOT::kSTLmultiset:
         case ROOT::kSTLmap:
         case ROOT::kSTLmultimap:
            if (narg>4 && strncmp(inside[4].c_str(),allocPrefix,allocPrefixLen)==0) {
               --narg;
            }
            break;
         case ROOT::kSTLunorderedset:
         case ROOT::kSTLunorderedmultiset:
            if (narg>5 && strncmp(inside[5].c_str(),allocPrefix,allocPrefixLen)==0) {
               --narg;
            }
            break;
         case ROOT::kSTLunorderedmap:
         case ROOT::kSTLunorderedmultimap:
            if (narg>6 && strncmp(inside[6].c_str(),allocPrefix,allocPrefixLen)==0) {
               --narg;
            }
            break;
      }
      if (stlkind!=0) {
         TClass *key = TClass::GetClass(inside[1].c_str());

         if (key) {
            // We only need to translate to a vector is the key is a class
            // (for which we do not know the sorting).
            std::string what;
            switch ( stlkind )  {
               case ROOT::kSTLmap:
               case ROOT::kSTLunorderedmap:
               case ROOT::kSTLmultimap:
               case ROOT::kSTLunorderedmultimap: {
                  what = "std::pair<";
                  what += inside[1];
                  what += ",";
                  what += inside[2];
                  if (what[what.size()-1]=='>') {
                     what += " >";
                  } else {
                     what += ">";
                  }
                  inside.clear();
                  inside.push_back("std::vector");
                  inside.push_back(what);
                  narg = 2;
                  break;
               }
               case ROOT::kSTLset:
               case ROOT::kSTLunorderedset:
               case ROOT::kSTLmultiset:
               case ROOT::kSTLunorderedmultiset:
                  inside[0] = "std::vector";
                  break;
            }
         }
         if (strncmp(inside[0].c_str(),"std::",5) != 0) {
            inside[0] = "std::" + inside[0];
         }
      } else {
         static const char *stlnames[] = { "pair", "greater", "less", "allocator" };
         for(unsigned int in = 0; in < sizeof(stlnames)/sizeof(stlnames[0]); ++in) {
            if (strncmp( inside[0].c_str(), stlnames[in], strlen(stlnames[in])) == 0 ) {
               inside[0] = "std::" + inside[0];
               break;
            }
         }
      }
      newname = inside[0];
      newname.Append("<");
      newname.Append(inside[1]);
      for(unsigned int j=2; j<narg; ++j) {
         if (!inside[j].empty()) {
            newname.Append(",");
            newname.Append(inside[j]);
         }
      }
      if (newname[newname.Length()-1]=='>') {
         newname.Append(" >");
      } else {
         newname.Append(">");
      }
      if (nestedLoc) newname.Append(inside[nestedLoc]);
   } else if ( newname == "string" ) {
      newname = "std::string";
   }
   return newname;
}
