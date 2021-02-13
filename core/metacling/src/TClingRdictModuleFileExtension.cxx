/// \file TClingRdictModuleFileExtension.cxx
///
/// \brief The file contains facilities to work with C++ module files extensions
///        used to store rdict files.
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date May, 2019
///
/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TClingRdictModuleFileExtension.h"

#include "TClingUtils.h"

#include "clang/Frontend/FrontendDiagnostic.h"
#include "clang/Lex/HeaderSearchOptions.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Sema/Sema.h"
#include "clang/Serialization/ASTReader.h"
#include "clang/Serialization/Module.h"

#include "llvm/ADT/Hashing.h"
#include "llvm/Bitcode/BitstreamWriter.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/Path.h"
#include "llvm/Support/raw_ostream.h"

#include <fstream>
#include <sstream>

/// Rdict module extension block name.
const std::string ROOT_CLING_RDICT_BLOCK_NAME = "root.cling.rdict";

/// Rdict module extension major version number.
constexpr uint16_t ROOT_CLING_RDICT_VERSION_MAJOR = 1;

/// Rdict module extension minor version number.
///
/// When the format changes IN ANY WAY, this number should be incremented.
constexpr uint16_t ROOT_CLING_RDICT_VERSION_MINOR = 1;


TClingRdictModuleFileExtension::Writer::~Writer() {}

void TClingRdictModuleFileExtension::Writer::writeExtensionContents(clang::Sema &SemaRef, llvm::BitstreamWriter &Stream)
{

   const clang::LangOptions &Opts = SemaRef.getLangOpts();
   const clang::Preprocessor &PP = SemaRef.getPreprocessor();

   llvm::StringRef CachePath = PP.getHeaderSearchInfo().getHeaderSearchOpts().ModuleCachePath;
   std::string RdictsStart = "lib" + Opts.CurrentModule + "_";
   const std::string RdictsEnd = "_rdict.pcm";

   using namespace llvm;
   using namespace clang;
   using namespace clang::serialization;
   // Write an abbreviation for this record.
   auto Abv = std::make_shared<BitCodeAbbrev>();
   Abv->Add(BitCodeAbbrevOp(FIRST_EXTENSION_RECORD_ID));
   Abv->Add(BitCodeAbbrevOp(BitCodeAbbrevOp::Blob));
   auto Abbrev = Stream.EmitAbbrev(std::move(Abv));
   auto Abv1 = std::make_shared<BitCodeAbbrev>();
   Abv1->Add(BitCodeAbbrevOp(FIRST_EXTENSION_RECORD_ID + 1));
   Abv1->Add(BitCodeAbbrevOp(BitCodeAbbrevOp::Blob));
   auto Abbrev1 = Stream.EmitAbbrev(std::move(Abv1));

   // Write a dict files into the extension block.
   std::error_code EC;
   for (llvm::sys::fs::directory_iterator DirIt(CachePath, EC), DirEnd; DirIt != DirEnd && !EC; DirIt.increment(EC)) {
      StringRef FilePath(DirIt->path());
      if (llvm::sys::fs::is_directory(FilePath))
         continue;
      StringRef FileName = llvm::sys::path::filename(FilePath);
      if (FileName.startswith(RdictsStart) && FileName.endswith(RdictsEnd)) {

         uint64_t Record[] = {FIRST_EXTENSION_RECORD_ID};
         Stream.EmitRecordWithBlob(Abbrev, Record, FileName);

         uint64_t Record1[] = {FIRST_EXTENSION_RECORD_ID + 1};
         std::ifstream fp(FilePath, std::ios::binary);
         std::ostringstream os;
         os << fp.rdbuf();
         Stream.EmitRecordWithBlob(Abbrev1, Record1, StringRef(os.str()));
         fp.close();

         EC = llvm::sys::fs::remove(FilePath);
         assert(!EC && "Unable to close _rdict file");
      }
   }
}

extern "C" void TCling__RegisterRdictForLoadPCM(const std::string &pcmFileNameFullPath, llvm::StringRef *pcmContent);

TClingRdictModuleFileExtension::Reader::Reader(clang::ModuleFileExtension *Ext, clang::ASTReader &Reader,
                                               clang::serialization::ModuleFile &Mod,
                                               const llvm::BitstreamCursor &InStream)
   : ModuleFileExtensionReader(Ext), Stream(InStream)
{
   // Read the extension block.
   llvm::SmallVector<uint64_t, 4> Record;
   llvm::StringRef CurrentRdictName;
   while (true) {
      llvm::BitstreamEntry Entry = Stream.advanceSkippingSubblocks();
      switch (Entry.Kind) {
      case llvm::BitstreamEntry::SubBlock:
      case llvm::BitstreamEntry::EndBlock:
      case llvm::BitstreamEntry::Error: return;

      case llvm::BitstreamEntry::Record: break;
      }

      Record.clear();
      llvm::StringRef Blob;
      unsigned RecCode = Stream.readRecord(Entry.ID, Record, &Blob);
      using namespace clang::serialization;
      switch (RecCode) {
      case FIRST_EXTENSION_RECORD_ID: {
         CurrentRdictName = Blob;
         break;
      }
      case FIRST_EXTENSION_RECORD_ID + 1: {
         // FIXME: Remove the string copy in fPendingRdicts.
         std::string ResolvedFileName
            = ROOT::TMetaUtils::GetRealPath(Mod.FileName);
         llvm::StringRef ModDir = llvm::sys::path::parent_path(ResolvedFileName);
         llvm::SmallString<255> FullRdictName = ModDir;
         llvm::sys::path::append(FullRdictName, CurrentRdictName);
         TCling__RegisterRdictForLoadPCM(FullRdictName.str(), &Blob);
         break;
      }
      }
   }
}

TClingRdictModuleFileExtension::Reader::~Reader() {}

TClingRdictModuleFileExtension::~TClingRdictModuleFileExtension() {}

clang::ModuleFileExtensionMetadata TClingRdictModuleFileExtension::getExtensionMetadata() const
{
   const std::string UserInfo = "";
   return {ROOT_CLING_RDICT_BLOCK_NAME, ROOT_CLING_RDICT_VERSION_MAJOR, ROOT_CLING_RDICT_VERSION_MINOR, UserInfo};
}

llvm::hash_code TClingRdictModuleFileExtension::hashExtension(llvm::hash_code Code) const
{
   Code = llvm::hash_combine(Code, ROOT_CLING_RDICT_BLOCK_NAME);
   Code = llvm::hash_combine(Code, ROOT_CLING_RDICT_VERSION_MAJOR);
   Code = llvm::hash_combine(Code, ROOT_CLING_RDICT_VERSION_MINOR);

   return Code;
}

std::unique_ptr<clang::ModuleFileExtensionWriter>
TClingRdictModuleFileExtension::createExtensionWriter(clang::ASTWriter &)
{
   return std::unique_ptr<clang::ModuleFileExtensionWriter>(new Writer(this));
}

std::unique_ptr<clang::ModuleFileExtensionReader>
TClingRdictModuleFileExtension::createExtensionReader(const clang::ModuleFileExtensionMetadata &Metadata,
                                                      clang::ASTReader &Reader, clang::serialization::ModuleFile &Mod,
                                                      const llvm::BitstreamCursor &Stream)
{
   return std::unique_ptr<clang::ModuleFileExtensionReader>(
      new TClingRdictModuleFileExtension::Reader(this, Reader, Mod, Stream));
}
