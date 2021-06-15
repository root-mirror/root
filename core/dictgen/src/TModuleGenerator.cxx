// @(#)root/utils:$Id$
// Author: Axel Naumann, 2-13-07-02
// Author: Danilo Piparo, 2013, 2014

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//
// PCM writer.
//
////////////////////////////////////////////////////////////////////////////////

#include "TModuleGenerator.h"

#include "TClingUtils.h"
#include "RConfigure.h"
#include <ROOT/RConfig.hxx>

#include "cling/Interpreter/CIFactory.h"
#include "clang/Basic/SourceManager.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Lex/HeaderSearch.h"
#include "clang/Lex/Preprocessor.h"
#include "llvm/Support/Path.h"

#include <map>

#ifndef R__WIN32
#include <unistd.h>
#endif

#include <iostream>

using namespace ROOT;
using namespace clang;

////////////////////////////////////////////////////////////////////////////////

TModuleGenerator::TModuleGenerator(CompilerInstance *CI,
                                   bool inlineInputHeaders,
                                   const std::string &shLibFileName,
                                   bool writeEmptyRootPCM):
   fCI(CI),
   fIsPCH(shLibFileName == "allDict.cxx"),
   fIsInPCH(writeEmptyRootPCM),
   fInlineInputHeaders(inlineInputHeaders),
   fDictionaryName(llvm::sys::path::stem(shLibFileName)),
   fDemangledDictionaryName(llvm::sys::path::stem(shLibFileName)),
   fModuleDirName(llvm::sys::path::parent_path(shLibFileName)),
   fErrorCount(0)
{
   // Need to resolve _where_ to create the pcm
   // We default in the lib subdirectory
   // otherwise we put it in the same directory as the dictionary file (for ACLiC)
   if (fModuleDirName.empty()) {
      fModuleDirName = "./";
   } else {
      fModuleDirName += "/";
   }

   fModuleFileName = fModuleDirName
                     + ROOT::TMetaUtils::GetModuleFileName(fDictionaryName.c_str());

   // Clean the dictionary name from characters which are not accepted in C++
   std::string tmpName = fDictionaryName;
   fDictionaryName.clear();
   ROOT::TMetaUtils::GetCppName(fDictionaryName, tmpName.c_str());

   // .pcm -> .pch
   if (IsPCH()) fModuleFileName[fModuleFileName.length() - 1] = 'h';

   // Add a random string to the filename to avoid races.
   auto makeTempFile = [&](const char *suffix) {
      llvm::SmallString<64> resultPath;
      std::string pattern = fModuleDirName + fDictionaryName + "%%%%%%%%%%" + suffix;
      llvm::sys::fs::createUniqueFile(pattern, resultPath); // NOTE: this creates the (empty) file
      // Return the full buffer, so caller can use `.c_str()` on the temporary.
      return resultPath;
   };

   fUmbrellaName = makeTempFile("_dictUmbrella.h").c_str();
   fContentName = makeTempFile("_dictContent.h").c_str();
}

TModuleGenerator::~TModuleGenerator()
{
   unlink(fUmbrellaName.c_str());
   unlink(fContentName.c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// Check whether the file's extension is compatible with C or C++.
/// Return whether source, header, Linkdef or nothing.

TModuleGenerator::ESourceFileKind
TModuleGenerator::GetSourceFileKind(const char *filename) const
{
   if (filename[0] == '-') return kSFKNotC;

   if (ROOT::TMetaUtils::IsLinkdefFile(filename)) {
      return kSFKLinkdef;
   }

   const size_t len = strlen(filename);
   const char *ext = filename + len - 1;
   while (ext >= filename && *ext != '.') --ext;
   if (ext < filename || *ext != '.') {
      // This might still be a system header, let's double check
      // via the FileManager.
      clang::Preprocessor &PP = fCI->getPreprocessor();
      clang::HeaderSearch &HdrSearch = PP.getHeaderSearchInfo();
      const clang::DirectoryLookup *CurDir = nullptr;
      const clang::FileEntry *hdrFileEntry
         =  HdrSearch.LookupFile(filename, clang::SourceLocation(),
                                 true /*isAngled*/, nullptr /*FromDir*/, CurDir,
                                 clang::ArrayRef<std::pair<const clang::FileEntry*,
                                                           const clang::DirectoryEntry*>>(),
                                 nullptr /*SearchPath*/,/*RelativePath*/ nullptr,
                                 nullptr /*RequestingModule*/, nullptr /*SuggestedModule*/,
                                 nullptr /*IsMapped*/, nullptr /*IsFrameworkFound*/);
      if (hdrFileEntry) {
         return kSFKHeader;
      }
      return kSFKNotC;
   }
   ++ext;
   const size_t lenExt = filename + len - ext;

   ESourceFileKind ret = kSFKNotC;
   switch (lenExt) {
      case 1: {
            const char last = toupper(filename[len - 1]);
            if (last == 'H') ret = kSFKHeader;
            else if (last == 'C') ret = kSFKSource;
            break;
         }
      case 2: {
            if (filename[len - 2] == 'h' && filename[len - 1] == 'h')
               ret = kSFKHeader;
            else if (filename[len - 2] == 'c' && filename[len - 1] == 'c')
               ret = kSFKSource;
            break;
         }
      case 3: {
            const char last = filename[len - 1];
            if ((last == 'x' || last == 'p')
                  && filename[len - 2] == last) {
               if (filename[len - 3] == 'h') ret = kSFKHeader;
               else if (filename[len - 3] == 'c') ret = kSFKSource;
            }
         }
   } // switch extension length

   return ret;
}


static
std::pair<std::string, std::string> SplitPPDefine(const std::string &in)
{
   std::string::size_type posEq = in.find('=');
   // No equal found: define to 1
   if (posEq == std::string::npos)
      return std::make_pair(in, "1");

   // Equal found
   return std::pair<std::string, std::string>
          (in.substr(0, posEq), in.substr(posEq + 1, std::string::npos));
}

////////////////////////////////////////////////////////////////////////////////
/// Parse -I -D -U headers.h SomethingLinkdef.h.

void TModuleGenerator::ParseArgs(const std::vector<std::string> &args)
{
   for (size_t iPcmArg = 1 /*skip argv0*/, nPcmArg = args.size();
         iPcmArg < nPcmArg; ++iPcmArg) {
      ESourceFileKind sfk = GetSourceFileKind(args[iPcmArg].c_str());
      if (sfk == kSFKHeader || sfk == kSFKSource) {
         fHeaders.push_back(args[iPcmArg]);
      } else if (sfk == kSFKLinkdef) {
         fLinkDefFile = args[iPcmArg];
      } else if (sfk == kSFKNotC && args[iPcmArg][0] == '-') {
         switch (args[iPcmArg][1]) {
            case 'I':
               if (args[iPcmArg] != "-I." &&  args[iPcmArg] != "-Iinclude") {
                  fCompI.push_back(args[iPcmArg].c_str() + 2);
               }
               break;
            case 'D':
               if (args[iPcmArg] != "-DTRUE=1" && args[iPcmArg] != "-DFALSE=0"
                     && args[iPcmArg] != "-DG__NOCINTDLL") {
                  fCompD.push_back(SplitPPDefine(args[iPcmArg].c_str() + 2));
               }
               break;
            case 'U':
               fCompU.push_back(args[iPcmArg].c_str() + 2);
               break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write
/// #ifndef FOO
/// # define FOO=bar
/// #endif

std::ostream &TModuleGenerator::WritePPDefines(std::ostream &out) const
{
   for (auto const & strPair : fCompD) {
      std::string cppname(strPair.first);
      size_t pos = cppname.find('(');
      if (pos != std::string::npos) cppname.erase(pos);
      out << "#ifndef " << cppname << "\n"
          "  #define " << strPair.first;
      out << " " << strPair.second;
      out << "\n"
          "#endif\n";
   }
   out << std::endl;
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// Write
/// #ifdef FOO
/// # undef FOO
/// #endif

std::ostream &TModuleGenerator::WritePPUndefines(std::ostream &out) const
{
   for (auto const & undef : fCompU) {
      out << "#ifdef " << undef << "\n"
          "  #undef " << undef << "\n"
          "#endif\n";
   }
   out << std::endl;
   return out;
}

////////////////////////////////////////////////////////////////////////////////
/// To be replaced with proper pragma handlers.

int WarnIfPragmaOnceDetected(const std::string& fullHeaderPath,
                              const std::string& headerFileContent)
{
   std::istringstream headerFile(headerFileContent);
   std::string line;
   while(std::getline(headerFile,line)){
      llvm::StringRef lineRef (line);
      auto trimmedLineRef = lineRef.trim();
      if (trimmedLineRef.startswith("#pragma") &&
          (trimmedLineRef.endswith(" once") || trimmedLineRef.endswith("\tonce"))) {
         std::cerr << "Error: #pragma once directive detected in header file "
                  << fullHeaderPath
                  << " which was requested to be inlined.\n";
         return 1;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

int ExtractBufferContent(const std::string& fullHeaderPath, std::string& bufferContent)
{
   std::ifstream buffer(fullHeaderPath);
   bufferContent = std::string((std::istreambuf_iterator<char>(buffer)),
                                std::istreambuf_iterator<char>());

   return WarnIfPragmaOnceDetected(fullHeaderPath,bufferContent);
}

////////////////////////////////////////////////////////////////////////////////
/// Write
/// #include "header1.h"
/// #include "header2.h"
/// or, if inlining of headers is requested, dump the content of the files.

std::ostream &TModuleGenerator::WritePPIncludes(std::ostream &out) const
{
   std::string fullHeaderPath;
   for (auto const & incl : fHeaders) {
      if (fInlineInputHeaders){
         bool headerFound = FindHeader(incl,fullHeaderPath);
         if (!headerFound){
            ROOT::TMetaUtils::Error(nullptr, "Cannot find header %s: cannot inline it.\n", fullHeaderPath.c_str());
            continue;
         }

         std::string bufferContent;
         fErrorCount += ExtractBufferContent(fullHeaderPath, bufferContent);

         out << bufferContent << std::endl;
      } else {
         out << "#include \"" << incl << "\"\n";
      }
   }
   out << std::endl;
   return out;
}

////////////////////////////////////////////////////////////////////////////////

std::ostream &TModuleGenerator::WriteStringVec(const std::vector<std::string> &vec,
      std::ostream &out) const
{
   for (auto const & theStr : vec) {
      out << "\"" << theStr << "\",\n";
   }
   out << "nullptr" << std::endl;
   return out;
}

////////////////////////////////////////////////////////////////////////////////

std::ostream &TModuleGenerator::WriteStringPairVec(const StringPairVec_t &vec,
      std::ostream &out) const
{
   for (auto const & strPair : vec) {
      out << "\"" << strPair.first;
      if (!strPair.second.empty()) {
         out << "=";
         // Need to escape the embedded quotes.
         for (const char *c = strPair.second.c_str(); *c != '\0'; ++c) {
            if (*c == '"') {
               out << "\\\"";
            } else {
               out << *c;
            }
         }
      }
      out << "\",\n";
   }
   out << "nullptr" << std::endl;
   return out;
}


void TModuleGenerator::WriteRegistrationSourceImpl(std::ostream& out,
                                                   const std::string &dictName,
                                                   const std::string &demangledDictName,
                                                   const std::vector<std::string> &headerArray,
                                                   const std::vector<std::string> &includePathArray,
                                                   const std::string &fwdDeclStringRAW,
                                                   const std::string &fwdDeclnArgsToKeepString,
                                                   const std::string &payloadCodeWrapped,
                                                   const std::string &headersClassesMapString,
                                                   const std::string &extraIncludes,
                                                   bool hasCxxModule) const
{
   // Dictionary initialization code for loading the module
   out << "namespace {\n"
          "  void TriggerDictionaryInitialization_" << dictName << "_Impl() {\n"
          "    static const char* headers[] = {\n";
   WriteStringVec(headerArray, out) << "    };\n";
   out << "    static const char* includePaths[] = {\n";
   WriteStringVec(includePathArray, out)
       << "    };\n";

   out << "    static const char* fwdDeclCode = " << fwdDeclStringRAW << ";\n"
       << "    static const char* payloadCode = " << payloadCodeWrapped << ";\n";
   // classesHeaders may depen on payloadCode
   out << "    static const char* classesHeaders[] = {\n"
       << headersClassesMapString
       << "\n};\n";
   out << "    static bool isInitialized = false;\n"
          "    if (!isInitialized) {\n"
          "      TROOT::RegisterModule(\"" << demangledDictName << "\",\n"
          "        headers, includePaths, payloadCode, fwdDeclCode,\n"
          "        TriggerDictionaryInitialization_" << dictName << "_Impl, "
                     << fwdDeclnArgsToKeepString << ", classesHeaders, "
                     << (hasCxxModule ? "/*hasCxxModule*/true" : "/*hasCxxModule*/false")
                     << ");\n"
          "      isInitialized = true;\n"
          "    }\n"
          "  }\n"
          "  static struct DictInit {\n"
          "    DictInit() {\n"
          "      TriggerDictionaryInitialization_" << dictName << "_Impl();\n"
          "    }\n"
          "    ~DictInit() {\n"
          "        TROOT::UnRegisterModule(\"" << demangledDictName << "\","
                                           "TriggerDictionaryInitialization_" << dictName << "_Impl);\n"
          "    }\n"
          "  } __TheDictionaryInitializer;\n"
          "}\n"
          "void TriggerDictionaryInitialization_" << dictName << "() {\n"
          "  TriggerDictionaryInitialization_" << dictName << "_Impl();\n"
          "}\n";
}

////////////////////////////////////////////////////////////////////////////////

void TModuleGenerator::WriteRegistrationSource(std::ostream &out, const std::string &fwdDeclnArgsToKeepString,
                                               const std::string &headersClassesMapString,
                                               const std::string &fwdDeclString, const std::string &extraIncludes, bool hasCxxModule) const
{
   if (hasCxxModule) {
      std::string emptyStr = "\"\"";
      WriteRegistrationSourceImpl(out, GetDictionaryName(), GetDemangledDictionaryName(), {}, {},
                                  fwdDeclString, "{}",
                                  emptyStr, headersClassesMapString, "0",
                                  /*HasCxxModule*/ true);
      return;
   }

   std::string fwdDeclStringSanitized = fwdDeclString;
#ifdef R__WIN32
   // Visual sudio has a limitation of 2048 characters max in raw strings, so split
   // the potentially huge DICTFWDDCLS raw string into multiple smaller ones
   constexpr char from[] = "\n";
   constexpr char to[] = "\n)DICTFWDDCLS\"\nR\"DICTFWDDCLS(";
   size_t start_pos = 0;
   while ((start_pos = fwdDeclStringSanitized.find(from, start_pos)) != std::string::npos) {
      if (fwdDeclStringSanitized.find(from, start_pos + 1) == std::string::npos) // skip the last
         break;
      if ((fwdDeclStringSanitized.at(start_pos + 1) == '}') || (fwdDeclStringSanitized.at(start_pos + 1) == '\n'))
         start_pos += 2;
      else {
         fwdDeclStringSanitized.replace(start_pos, strlen(from), to);
         start_pos += strlen(to); // In case 'to' contains 'from', like replacing 'x' with 'yx'
      }
   }
#endif
   std::string fwdDeclStringRAW;
   if ("nullptr" == fwdDeclStringSanitized || "\"\"" == fwdDeclStringSanitized) {
      fwdDeclStringRAW = fwdDeclStringSanitized;
   } else {
      fwdDeclStringRAW = "R\"DICTFWDDCLS(\n";
      fwdDeclStringRAW += "#line 1 \"";
      fwdDeclStringRAW += fDictionaryName +" dictionary forward declarations' payload\"\n";
      fwdDeclStringRAW += fwdDeclStringSanitized;
      fwdDeclStringRAW += ")DICTFWDDCLS\"";
   }

   if (fIsInPCH)
      fwdDeclStringRAW = "nullptr";

   std::string payloadCode;

   // Increase the value of the diagnostics pointing out from which
   // dictionary this payload comes from in case of errors
   payloadCode += "#line 1 \"";
   payloadCode += fDictionaryName +" dictionary payload\"\n";

   // Add defines and undefines to the payloadCode
   std::ostringstream definesAndUndefines;
   //    Anticipate the undefines.
   //    Suppose to have a namespace called "declarations" used in R5 for template
   //    instantiations in the header given to genreflex.
   //    Now, in this namespace, objects with some names, typically dummy, will be
   //    present.
   //    If you give such headers to cling to parse, problems will occour, as the
   //    names appear multiple times. One possible solution is to get out of this
   //    with preprocessor defines given to genreflex, redefining "declarations"
   //    to a hash or <project>_<package> via the build system.
   WritePPUndefines(definesAndUndefines);
   WritePPDefines(definesAndUndefines);
   payloadCode += definesAndUndefines.str();

   // If necessary, inline the headers
   std::string hdrFullPath;
   std::string inlinedHeaders;

   auto findAndAddToInlineHeaders = [&](const std::string& hdrName) {
      bool headerFound = FindHeader(hdrName,hdrFullPath);
      if (!headerFound) {
         ROOT::TMetaUtils::Error(nullptr, "Cannot find header %s: cannot inline it.\n", hdrName.c_str());
      } else {
         std::ifstream headerFile(hdrFullPath.c_str());
         const std::string headerFileAsStr((std::istreambuf_iterator<char>(headerFile)),
                                            std::istreambuf_iterator<char>());
         inlinedHeaders += headerFileAsStr;
      }
   };

   if (fInlineInputHeaders) {
      for (auto const & hdrName : fHeaders) {
        findAndAddToInlineHeaders(hdrName);
      }

   } else {
      // Now, if not, just #include them in the payload except for the linkdef
      // file, if any.
      for (auto & hdrName : fHeaders) {
         inlinedHeaders += "#include \"" + hdrName + "\"\n";
      }
   }

   if (0 != fLinkDefFile.size() && !fIsPCH) {
      findAndAddToInlineHeaders(fLinkDefFile);
   }

   // Recover old genreflex behaviour, i.e. do not print warnings due to glitches
   // in the headers at runtime. This is not synonym of ignoring warnings as they
   // will be printed at dictionary generation time.
   // In order to do this we leverage the diagnostic pragmas and, since there is no
   // way to express as a pragma the option "-Wno-deprecated" the
   // _BACKWARD_BACKWARD_WARNING_H macro, used to avoid to go through
   // backward/backward_warning.h.
   payloadCode += "#define _BACKWARD_BACKWARD_WARNING_H\n"
                  "// Inline headers\n"+
                  inlinedHeaders + "\n"+
                  (extraIncludes.empty() ? "" : "// Extra includes\n" + extraIncludes + "\n") +
                  "#undef  _BACKWARD_BACKWARD_WARNING_H\n";

   // We cannot stream the contents in strings and pass it to
   // WriteRegistrationSourceImpl because we exceed the std::string::max_size on
   // Windows.
   std::vector<std::string> headerArray = {"0"};
   if (!fInlineInputHeaders)
      headerArray = fHeaders;
   const std::vector<std::string>& includePathArray = fCompI;

   std::string payloadcodeWrapped = "nullptr";
   if (!fIsInPCH)
      payloadcodeWrapped = "R\"DICTPAYLOAD(\n" + payloadCode + ")DICTPAYLOAD\"";

   WriteRegistrationSourceImpl(out, GetDictionaryName(),
                               GetDemangledDictionaryName(),
                               headerArray,
                               includePathArray,
                               fwdDeclStringRAW,
                               fwdDeclnArgsToKeepString,
                               payloadcodeWrapped,
                               headersClassesMapString,
                               extraIncludes,
                               /*HasCxxModule*/false);
}

////////////////////////////////////////////////////////////////////////////////
/// Write a header file describing the content of this module
/// through a series of variables inside the namespace
/// ROOT::Dict::[DictionaryName]. Each variable is an array of string
/// literals, with a const char* of 0 being the last element, e.g.
/// ROOT::Dict::_DictName::arrIncludes[] = { "A.h", "B.h", 0 };

void TModuleGenerator::WriteContentHeader(std::ostream &out) const
{
   out << "namespace ROOT { namespace Dict { namespace _"
       << GetDictionaryName() << "{\n";

   out << "const char* arrIncludes[] = {\n";
   WriteHeaderArray(out) << "};\n";

   out << "const char* arrIncludePaths[] = {\n";
   WriteIncludePathArray(out) << "};\n";
   /*
      out << "const char* arrDefines[] = {\n";
      WriteDefinesArray(out) << "};\n";

      out << "const char* arrUndefines[] = {\n";
      WriteUndefinesArray(out) << "};\n";*/

   out << "} } }" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the header is found in the include paths
/// in this case also fill the full path variable with the full path.

bool TModuleGenerator::FindHeader(const std::string &hdrName, std::string &hdrFullPath) const
{
   hdrFullPath = hdrName;
   if (llvm::sys::fs::exists(hdrFullPath))
      return true;
   for (auto const &incDir : fCompI) {
      hdrFullPath = incDir + ROOT::TMetaUtils::GetPathSeparator() + hdrName;
      if (llvm::sys::fs::exists(hdrFullPath)) {
         return true;
      }
   }
   clang::Preprocessor &PP = fCI->getPreprocessor();
   clang::HeaderSearch &HdrSearch = PP.getHeaderSearchInfo();
   const clang::DirectoryLookup *CurDir = nullptr;
   if (const clang::FileEntry *hdrFileEntry
         =  HdrSearch.LookupFile(hdrName, clang::SourceLocation(),
                                 true /*isAngled*/, nullptr /*FromDir*/, CurDir,
                                 clang::ArrayRef<std::pair<const clang::FileEntry*,
                                                         const clang::DirectoryEntry*>>(),
                                 nullptr /*SearchPath*/, nullptr /*RelativePath*/,
                                 nullptr /*RequestingModule*/, nullptr/*SuggestedModule*/,
                                 nullptr /*IsMapped*/, nullptr /*IsFrameworkFound*/)) {
      hdrFullPath = hdrFileEntry->getName();
      return true;
   }

   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a header file pulling in the content of this module
/// through a series of #defined, #undefs and #includes.
/// The sequence corrsponds to a rootcling invocation with
///   -c -DFOO -UBAR header.h
/// I.e. defines, undefines and finally includes.

void TModuleGenerator::WriteUmbrellaHeader(std::ostream &out) const
{
   WritePPDefines(out);
   WritePPUndefines(out);
   WritePPIncludes(out);
}
