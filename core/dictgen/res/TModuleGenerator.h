// @(#)root/utils:$Id$
// Author: Axel Naumann, 2-13-07-02

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
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

#ifndef ROOT_TModuleGenerator
#define ROOT_TModuleGenerator

#include <ostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>

namespace clang {
   class CompilerInstance;
   class SourceManager;
}

namespace ROOT {

//______________________________________________________________________________
   class TModuleGenerator {
   public:
      enum ESourceFileKind {
         kSFKNotC,
         kSFKHeader,
         kSFKSource,
         kSFKLinkdef
      };

      TModuleGenerator(clang::CompilerInstance *CI,
                       bool inlineHeader,
                       const std::string &shLibFileName,
                       bool isInPCH);
      ~TModuleGenerator();

      // FIXME: remove once PCH is gone.
      bool IsPCH() const {
         return fIsPCH;
      }
      void ParseArgs(const std::vector<std::string> &args);

      const std::string &GetDictionaryName() const {
         return fDictionaryName;
      }

      const std::string &GetDemangledDictionaryName() const {
         return fDemangledDictionaryName;
      }

      const std::string &GetModuleFileName() const {
         return fModuleFileName;
      }
      const std::string &GetModuleDirName() const {
         return fModuleDirName;
      }

      int GetErrorCount() const {
         return fErrorCount;
      }

      const std::string &GetUmbrellaName() const {
         return fUmbrellaName;
      }
      const std::string &GetContentName() const {
         return fContentName;
      }

      const std::vector<std::string> &GetHeaders() const {
         return fHeaders;
      }

      const std::vector<std::string> &GetIncludePaths() const {
         return fCompI;
      }

      std::ostream &WritePPDefines(std::ostream &out) const;
      std::ostream &WritePPUndefines(std::ostream &out) const;

      void WriteRegistrationSource(std::ostream &out, const std::string &fwdDeclnArgsToKeepString,
                                   const std::string &headersClassesMapString, const std::string &fwdDeclsString,
                                   const std::string &extraIncludes, bool hasCxxModule) const;
      void WriteContentHeader(std::ostream &out) const;
      void WriteUmbrellaHeader(std::ostream &out) const;

   private:
      void WriteRegistrationSourceImpl(std::ostream& out,
                                       const std::string &dictName,
                                       const std::string &demangledDictName,
                                       const std::vector<std::string> &headerArray,
                                       const std::vector<std::string> &includePathArray,
                                       const std::string &fwdDeclStringRAW,
                                       const std::string &fwdDeclnArgsToKeepString,
                                       const std::string &payloadCodeWrapped,
                                       const std::string &headersClassesMapString,
                                       const std::string &extraIncludes,
                                       bool hasCxxModule) const;

      void ConvertToCppString(std::string &text) const;

      std::ostream &WritePPIncludes(std::ostream &out) const;

      std::ostream &WritePPCode(std::ostream &out) const {
         // Write defines, undefiles, includes, corrsponding to a rootcling
         // invocation with -c -DFOO -UBAR header.h.
         WritePPDefines(out);
         WritePPUndefines(out);
         return WritePPIncludes(out);
      }

      std::ostream &WriteHeaderArray(std::ostream &out) const {
         // Write "header1.h",\n"header2.h",\n0\n
         return WriteStringVec(fHeaders, out);
      }
      std::ostream &WriteIncludePathArray(std::ostream &out) const {
         // Write "./include",\n"/usr/include",\n0\n
         return WriteStringVec(fCompI, out);
      }
      std::ostream &WriteDefinesArray(std::ostream &out) const {
         // Write "DEFINED",\n"NAME=\"VALUE\"",\n0\n
         return WriteStringPairVec(fCompD, out);
      }
      std::ostream &WriteUndefinesArray(std::ostream &out) const {
         // Write "UNDEFINED",\n"NAME",\n0\n
         return WriteStringVec(fCompU, out);
      }

      bool FindHeader(const std::string& hdrName, std::string& hdrFullPath) const;

      typedef std::vector<std::pair<std::string, std::string> > StringPairVec_t;

      ESourceFileKind GetSourceFileKind(const char *filename) const;
      std::ostream &WriteStringVec(const std::vector<std::string> &vec,
                                   std::ostream &out) const;
      std::ostream &WriteStringPairVec(const StringPairVec_t &vecP,
                                       std::ostream &out) const;

      clang::CompilerInstance *fCI;
      bool fIsPCH;
      bool fIsInPCH; // whether the headers of this module are part of the PCH.
      bool fInlineInputHeaders;

      std::string fDictionaryName; // Name of the dictionary, e.g. "Base"
      std::string fDemangledDictionaryName; // Demangled name of the dictionary
      std::string fModuleFileName; // PCM file name
      std::string fModuleDirName; // PCM output directory
      std::string fUmbrellaName; // name of umbrella header in PCM
      std::string fContentName; // name of content description header in PCM

      std::vector<std::string> fHeaders; // exported headers in PCM
      std::string fLinkDefFile; // The name of the linkdef file
      std::vector<std::string> fCompI; // -I; needed only for ACLiC without PCMs

      StringPairVec_t fCompD; // -Dfirst=second
      std::vector<std::string> fCompU; // -Ufirst
      mutable int fErrorCount;
   };

} // namespace ROOT

#endif // ROOT_TModuleGenerator
