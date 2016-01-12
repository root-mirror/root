/// \file TFile.cxx
/// \ingroup Base ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-31
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/TFile.h"
#include "TFile.h"

#include <mutex>

ROOT::v7::TDirectory& ROOT::v7::TDirectory::Heap() {
  static TDirectory heapDir;
  return heapDir;
}

namespace {
/// We cannot afford users not closing their files. Yes, we return a unique_ptr -
/// but that might be stored in an object that itself leaks. That would leave
/// the TFile unclosed and data corrupted / not written. Instead, keep a
/// collection of all opened writable TFiles and close them at destruction time,
/// explicitly.
static void AddFilesToClose(ROOT::v7::TCoopPtr<ROOT::v7::Internal::TFileImplBase> pFile) {
  struct CloseFiles_t {
    std::vector<ROOT::v7::TCoopPtr<ROOT::v7::Internal::TFileImplBase>> fFiles;
    std::mutex fMutex;
    ~CloseFiles_t() {
      for (auto& pFile: fFiles)
        if (pFile)
          pFile->Flush(); // or Close()? but what if there's still a Write()?
    }
  };
  static CloseFiles_t closer;

  std::lock_guard<std::mutex> lock(closer.fMutex);
  closer.fFiles.emplace_back(pFile);
}

/** \class TFSFile
 TFileImplBase for a file-system (POSIX) style TFile.
 */
class TFileSystemFile: public ROOT::v7::Internal::TFileImplBase {
  ::TFile* fOldFile;

public:
  TFileSystemFile(const std::string& name, const char* mode):
    fOldFile(::TFile::Open(name.c_str(), mode)) {
  }

  void Flush() final { fOldFile->Flush(); }

  void Close() final { fOldFile->Close(); }

  ~TFileSystemFile() {
    delete fOldFile;
  }
};
}

ROOT::v7::TFilePtr::TFilePtr(ROOT::v7::TCoopPtr<ROOT::v7::Internal::TFileImplBase> impl):
fImpl(impl)
{
  AddFilesToClose(impl);
}


ROOT::v7::TFilePtr ROOT::v7::TFilePtr::OpenForRead(std::string_view name) {
  // will become delegation to TFileSystemFile, TWebFile etc.
  return TFilePtr(MakeCoop<TFileSystemFile>(name.to_string(), "READ"));
}
ROOT::v7::TFilePtr ROOT::v7::TFilePtr::Create(std::string_view name) {
  // will become delegation to TFileSystemFile, TWebFile etc.
  return TFilePtr(MakeCoop<TFileSystemFile>(name.to_string(), "CREATE"));
}
ROOT::v7::TFilePtr ROOT::v7::TFilePtr::Recreate(std::string_view name) {
  // will become delegation to TFileSystemFile, TWebFile etc.
  return TFilePtr(MakeCoop<TFileSystemFile>(name.to_string(), "RECREATE"));
}
ROOT::v7::TFilePtr ROOT::v7::TFilePtr::OpenForUpdate(std::string_view name) {
  // will become delegation to TFileSystemFile, TWebFile etc.
  return TFilePtr(MakeCoop<TFileSystemFile>(name.to_string(), "UPDATE"));
}
