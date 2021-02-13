// @(#)root/io:$Id$
// Author: Andreas Peters + Fons Rademakers + Rene Brun  26/5/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\class TFileMerger TFileMerger.cxx
\ingroup IO

This class provides file copy and merging services.

It can be used to copy files (not only ROOT files), using TFile or
any of its remote file access plugins. It is therefore useful in
a Grid environment where the files might be accessible only remotely.
The merging interface allows files containing histograms and trees
to be merged, like the standalone hadd program.
*/

#include "TFileMerger.h"
#include "TDirectory.h"
#include "TUrl.h"
#include "TFile.h"
#include "TUUID.h"
#include "TSystem.h"
#include "TKey.h"
#include "THashList.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TFileMergeInfo.h"
#include "TClassRef.h"
#include "TROOT.h"
#include "TMemFile.h"
#include "TVirtualMutex.h"

#ifdef WIN32
// For _getmaxstdio
#include <cstdio>
#else
// For getrlimit
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <cstring>

ClassImp(TFileMerger);

TClassRef R__TH1_Class("TH1");
TClassRef R__TTree_Class("TTree");

static const Int_t kCpProgress = BIT(14);
static const Int_t kCintFileNumber = 100;
////////////////////////////////////////////////////////////////////////////////
/// Return the maximum number of allowed opened files minus some wiggle room
/// for CINT or at least of the standard library (stdio).

static Int_t R__GetSystemMaxOpenedFiles()
{
   int maxfiles;
#ifdef WIN32
   maxfiles = _getmaxstdio();
#else
   rlimit filelimit;
   if (getrlimit(RLIMIT_NOFILE,&filelimit)==0) {
      maxfiles = filelimit.rlim_cur;
   } else {
      // We could not get the value from getrlimit, let's return a reasonable default.
      maxfiles = 512;
   }
#endif
   if (maxfiles > kCintFileNumber) {
      return maxfiles - kCintFileNumber;
   } else if (maxfiles > 5) {
      return maxfiles - 5;
   } else {
      return maxfiles;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create file merger object.

TFileMerger::TFileMerger(Bool_t isLocal, Bool_t histoOneGo)
            : fMaxOpenedFiles( R__GetSystemMaxOpenedFiles() ),
              fLocal(isLocal), fHistoOneGo(histoOneGo)
{
   fMergeList.SetOwner(kTRUE);
   fExcessFiles.SetOwner(kTRUE);

   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfCleanups()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TFileMerger::~TFileMerger()
{
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfCleanups()->Remove(this);
   }
   SafeDelete(fOutputFile);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset merger file list.

void TFileMerger::Reset()
{
   fFileList.Clear();
   fMergeList.Clear();
   fExcessFiles.Clear();
   fObjectNames.Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Add file to file merger.

Bool_t TFileMerger::AddFile(const char *url, Bool_t cpProgress)
{
   if (fPrintLevel > 0) {
      Printf("%s Source file %d: %s", fMsgPrefix.Data(), fFileList.GetEntries() + fExcessFiles.GetEntries() + 1, url);
   }

   TFile *newfile = 0;
   TString localcopy;

   if (fFileList.GetEntries() >= (fMaxOpenedFiles-1)) {

      TObjString *urlObj = new TObjString(url);
      fMergeList.Add(urlObj);

      urlObj = new TObjString(url);
      urlObj->SetBit(kCpProgress);
      fExcessFiles.Add(urlObj);
      return kTRUE;
   }

   // We want gDirectory untouched by anything going on here
   TDirectory::TContext ctxt;

   if (fLocal) {
      TUUID uuid;
      localcopy.Form("file:%s/ROOTMERGE-%s.root", gSystem->TempDirectory(), uuid.AsString());
      if (!TFile::Cp(url, localcopy, cpProgress)) {
         Error("AddFile", "cannot get a local copy of file %s", url);
         return kFALSE;
      }
      newfile = TFile::Open(localcopy, "READ");
   } else {
      newfile = TFile::Open(url, "READ");
   }

   // Zombie files should also be skipped
   if (newfile && newfile->IsZombie()) {
      delete newfile;
      newfile = 0;
   }

   if (!newfile) {
      if (fLocal)
         Error("AddFile", "cannot open local copy %s of URL %s",
               localcopy.Data(), url);
      else
         Error("AddFile", "cannot open file %s", url);
      return kFALSE;
   } else {
      if (fOutputFile && fOutputFile->GetCompressionLevel() != newfile->GetCompressionLevel()) fCompressionChange = kTRUE;

      newfile->SetBit(kCanDelete);
      fFileList.Add(newfile);

      TObjString *urlObj = new TObjString(url);
      fMergeList.Add(urlObj);

      return  kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add the TFile to this file merger and *do not* give ownership of the TFile to this
/// object.
///
/// Return kTRUE if the addition was successful.

Bool_t TFileMerger::AddFile(TFile *source, Bool_t cpProgress)
{
   return AddFile(source,kFALSE,cpProgress);
}

////////////////////////////////////////////////////////////////////////////////
/// Add the TFile to this file merger and give ownership of the TFile to this
/// object (unless kFALSE is returned).
///
/// Return kTRUE if the addition was successful.

Bool_t TFileMerger::AddAdoptFile(TFile *source, Bool_t cpProgress)
{
   return AddFile(source,kTRUE,cpProgress);
}

////////////////////////////////////////////////////////////////////////////////
/// Add the TFile to this file merger and give ownership of the TFile to this
/// object (unless kFALSE is returned).
///
/// Return kTRUE if the addition was successful.

Bool_t TFileMerger::AddFile(TFile *source, Bool_t own, Bool_t cpProgress)
{
   if (source == 0 || source->IsZombie()) {
      return kFALSE;
   }

   if (fPrintLevel > 0) {
      Printf("%s Source file %d: %s",fMsgPrefix.Data(),fFileList.GetEntries()+1,source->GetName());
   }

   TFile *newfile = 0;
   TString localcopy;

   // We want gDirectory untouched by anything going on here
   TDirectory::TContext ctxt;
   if (fLocal && !source->InheritsFrom(TMemFile::Class())) {
      TUUID uuid;
      localcopy.Form("file:%s/ROOTMERGE-%s.root", gSystem->TempDirectory(), uuid.AsString());
      if (!source->Cp(localcopy, cpProgress)) {
         Error("AddFile", "cannot get a local copy of file %s", source->GetName());
         return kFALSE;
      }
      newfile = TFile::Open(localcopy, "READ");
      // Zombie files should also be skipped
      if (newfile && newfile->IsZombie()) {
         delete newfile;
         newfile = 0;
      }
   } else {
      newfile = source;
   }

   if (!newfile) {
      if (fLocal)
         Error("AddFile", "cannot open local copy %s of URL %s",
               localcopy.Data(), source->GetName());
      else
         Error("AddFile", "cannot open file %s", source->GetName());
      return kFALSE;
   } else {
      if (fOutputFile && fOutputFile->GetCompressionSettings() != newfile->GetCompressionSettings()) fCompressionChange = kTRUE;

      if (own || newfile != source) {
         newfile->SetBit(kCanDelete);
      } else {
         newfile->ResetBit(kCanDelete);
      }
      fFileList.Add(newfile);

      TObjString *urlObj = new TObjString(source->GetName());
      fMergeList.Add(urlObj);

      if (newfile != source && own) {
         delete source;
      }
      return  kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Open merger output file.

Bool_t TFileMerger::OutputFile(const char *outputfile, Bool_t force, Int_t compressionLevel)
{
   return OutputFile(outputfile,(force?"RECREATE":"CREATE"),compressionLevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Open merger output file.

Bool_t TFileMerger::OutputFile(const char *outputfile, Bool_t force)
{
   Bool_t res = OutputFile(outputfile,(force?"RECREATE":"CREATE"),1); // 1 is the same as the default from the TFile constructor.
   fExplicitCompLevel = kFALSE;
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Open merger output file.
///
/// The 'mode' parameter is passed to the TFile constructor as the option, it
/// should be one of 'NEW','CREATE','RECREATE','UPDATE'
/// 'UPDATE' is usually used in conjunction with IncrementalMerge.

Bool_t TFileMerger::OutputFile(const char *outputfile, const char *mode, Int_t compressionLevel)
{
   // We want gDirectory untouched by anything going on here
   TDirectory::TContext ctxt;
   if (TFile *outputFile = TFile::Open(outputfile, mode, "", compressionLevel))
      return OutputFile(std::unique_ptr<TFile>(outputFile));

   Error("OutputFile", "cannot open the MERGER output file %s", fOutputFilename.Data());
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set an output file opened externally by the users

Bool_t TFileMerger::OutputFile(std::unique_ptr<TFile> outputfile)
{
   if (!outputfile || outputfile->IsZombie()) {
      Error("OutputFile", "cannot open the MERGER output file %s", (outputfile) ? outputfile->GetName() : "");
      return kFALSE;
   }

   if (!outputfile->IsWritable()) {
      Error("OutputFile", "output file %s is not writable", outputfile->GetName());
      return kFALSE;
   }

   fExplicitCompLevel = kTRUE;

   TFile *oldfile = fOutputFile;
   fOutputFile = 0; // This avoids the complaint from RecursiveRemove about the file being deleted which is here
                    // spurrious. (see RecursiveRemove).
   SafeDelete(oldfile);

   fOutputFilename = outputfile->GetName();
   // We want gDirectory untouched by anything going on here
   TDirectory::TContext ctxt;
   fOutputFile = outputfile.release(); // Transfer the ownership of the file.

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Open merger output file.  'mode' is passed to the TFile constructor as the option, it should
/// be one of 'NEW','CREATE','RECREATE','UPDATE'
/// 'UPDATE' is usually used in conjunction with IncrementalMerge.

Bool_t TFileMerger::OutputFile(const char *outputfile, const char *mode /* = "RECREATE" */)
{
   Bool_t res = OutputFile(outputfile,mode,1); // 1 is the same as the default from the TFile constructor.
   fExplicitCompLevel = kFALSE;
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Print list of files being merged.

void TFileMerger::PrintFiles(Option_t *options)
{
   fFileList.Print(options);
   fExcessFiles.Print(options);
}

////////////////////////////////////////////////////////////////////////////////
/// Merge the files.
///
/// If no output file was specified it will write into
/// the file "FileMerger.root" in the working directory. Returns true
/// on success, false in case of error.

Bool_t TFileMerger::Merge(Bool_t)
{
   return PartialMerge(kAll | kRegular);
}

namespace {

/// Merge a list of RNTuples
Long64_t MergeRNTuples(TClass* rntupleHandle, const TString& /* target */, const TList& /* sources */) {
   if (!rntupleHandle) {
      return Long64_t(-1);
   }
   // todo(max) implement rntuple merger
   // [ ] build complete list of sources (some sources may actually be a directory with RNTuples inside)
   // [ ] merge them
   ROOT::MergeFunc_t func = rntupleHandle->GetMerge();
   return func(static_cast<void*>(rntupleHandle), nullptr, nullptr);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
/// Merge all objects in a directory
///
/// The type is defined by the bit values in TFileMerger::EPartialMergeType.

Bool_t TFileMerger::MergeRecursive(TDirectory *target, TList *sourcelist, Int_t type /* = kRegular | kAll */)
{
   Bool_t status = kTRUE;
   Bool_t onlyListed = kFALSE;
   if (fPrintLevel > 0) {
      Printf("%s Target path: %s",fMsgPrefix.Data(),target->GetPath());
   }

   // Get the dir name
   TString path(target->GetPath());
   // coverity[unchecked_value] 'target' is from a file so GetPath always returns path starting with filename:
   path.Remove(0, std::strlen(target->GetFile()->GetPath()));

   Int_t nguess = sourcelist->GetSize()+1000;
   THashList allNames(nguess);
   allNames.SetOwner(kTRUE);
   // If the mode is set to skipping list objects, add names to the allNames list
   if (type & kSkipListed) {
      TObjArray *arr = fObjectNames.Tokenize(" ");
      arr->SetOwner(kFALSE);
      for (Int_t iname=0; iname<arr->GetEntriesFast(); iname++)
         allNames.Add(arr->At(iname));
      delete arr;
   }
   ((THashList*)target->GetList())->Rehash(nguess);
   ((THashList*)target->GetListOfKeys())->Rehash(nguess);

   TFileMergeInfo info(target);
   info.fIOFeatures = fIOFeatures;
   info.fOptions = fMergeOptions;
   if (fFastMethod && ((type&kKeepCompression) || !fCompressionChange) ) {
      info.fOptions.Append(" fast");
   }

   TFile      *current_file;
   TDirectory *current_sourcedir;
   if (type & kIncremental) {
      current_file      = 0;
      current_sourcedir = target;
   } else {
      current_file      = (TFile*)sourcelist->First();
      current_sourcedir = current_file->GetDirectory(path);
   }
   while (current_file || current_sourcedir) {
      // When current_sourcedir != 0 and current_file == 0 we are going over the target
      // for an incremental merge.
      if (current_sourcedir && (current_file == 0 || current_sourcedir != target)) {

         // loop over all keys in this directory
         TIter nextkey( current_sourcedir->GetListOfKeys() );
         TKey *key;
         TString oldkeyname;

         while ( (key = (TKey*)nextkey())) {

            // Keep only the highest cycle number for each key for mergeable objects. They are stored
            // in the (hash) list consecutively and in decreasing order of cycles, so we can continue
            // until the name changes. We flag the case here and we act consequently later.
            Bool_t alreadyseen = (oldkeyname == key->GetName()) ? kTRUE : kFALSE;

            // Read in but do not copy directly the processIds.
            if (strcmp(key->GetClassName(),"TProcessID") == 0) { key->ReadObj(); continue;}

            // If we have already seen this object [name], we already processed
            // the whole list of files for this objects and we can just skip it
            // and any related cycles.
            if (allNames.FindObject(key->GetName())) {
               oldkeyname = key->GetName();
               continue;
            }

            TClass *cl = TClass::GetClass(key->GetClassName());
            if (!cl) {
               Info("MergeRecursive", "cannot indentify object type (%s), name: %s title: %s",
                    key->GetClassName(), key->GetName(), key->GetTitle());
               continue;
            }
            // For mergeable objects we add the names in a local hashlist handling them
            // again (see above)
            if (cl->GetMerge() || cl->InheritsFrom(TDirectory::Class()) ||
               (cl->IsTObject() && !cl->IsLoaded() &&
                 /* If it has a dictionary and GetMerge() is nullptr then we already know the answer
                    to the next question is 'no, if we were to ask we would useless trigger
                    auto-parsing */
                 (cl->GetMethodWithPrototype("Merge", "TCollection*,TFileMergeInfo*") ||
                  cl->GetMethodWithPrototype("Merge", "TCollection*"))))
               allNames.Add(new TObjString(key->GetName()));

            if (fNoTrees && cl->InheritsFrom(R__TTree_Class)) {
               // Skip the TTree objects and any related cycles.
               oldkeyname = key->GetName();
               continue;
            }
            // Check if only the listed objects are to be merged
            if (type & kOnlyListed) {
               onlyListed = kFALSE;
               oldkeyname = key->GetName();
               oldkeyname += " ";
               onlyListed = fObjectNames.Contains(oldkeyname);
               oldkeyname = key->GetName();
               if ((!onlyListed) && (!cl->InheritsFrom(TDirectory::Class()))) continue;
            }

            if (!(type&kResetable && type&kNonResetable)) {
               // If neither or both are requested at the same time, we merger both types.
               if (!(type&kResetable)) {
                  if (cl->GetResetAfterMerge()) {
                     // Skip the object with a reset after merge routine (TTree and other incrementally mergeable objects)
                     oldkeyname = key->GetName();
                     continue;
                  }
               }
               if (!(type&kNonResetable)) {
                  if (!cl->GetResetAfterMerge()) {
                     // Skip the object without a reset after merge routine (Histograms and other non incrementally mergeable objects)
                     oldkeyname = key->GetName();
                     continue;
                  }
               }
            }
            // read object from first source file
            TObject *obj;
            if (type & kIncremental) {
               obj = current_sourcedir->GetList()->FindObject(key->GetName());
               if (!obj) {
                  obj = key->ReadObj();
               }
            } else {
               obj = key->ReadObj();
            }
            if (!obj) {
               Info("MergeRecursive", "could not read object for key {%s, %s}",
                    key->GetName(), key->GetTitle());
               continue;
            }
            // if (cl->IsTObject())
            //    obj->ResetBit(kMustCleanup);
            if (cl->IsTObject() && cl != obj->IsA()) {
               Error("MergeRecursive", "TKey and object retrieve disagree on type (%s vs %s).  Continuing with %s.",
                    key->GetClassName(), obj->IsA()->GetName(), obj->IsA()->GetName());
               cl = obj->IsA();
            }
            Bool_t canBeMerged = kTRUE;

            if ( cl->InheritsFrom( TDirectory::Class() ) ) {
               // it's a subdirectory

               target->cd();
               TDirectory *newdir;

               // For incremental or already seen we may have already a directory created
               if (type & kIncremental || alreadyseen) {
                  newdir = target->GetDirectory(obj->GetName());
                  if (!newdir) {
                     newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
                     // newdir->ResetBit(kMustCleanup);
                  }
               } else {
                  newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
                  // newdir->ResetBit(kMustCleanup);
               }

               // newdir is now the starting point of another round of merging
               // newdir still knows its depth within the target file via
               // GetPath(), so we can still figure out where we are in the recursion

               // If this folder is a onlyListed object, merge everything inside.
               if (onlyListed) type &= ~kOnlyListed;
               status = MergeRecursive(newdir, sourcelist, type);
               if (onlyListed) type |= kOnlyListed;
               if (!status) return status;
            } else if (!cl->IsTObject() && cl->GetMerge()) {
               // merge objects that don't derive from TObject
               if (std::string(key->GetClassName()) == "ROOT::Experimental::RNTuple") {
                  Warning("MergeRecursive", "merging RNTuples is experimental");
                  // todo(max): check if this works when a TDirectory is passed as the first
                  // input argument
                  Long64_t mergeResult = MergeRNTuples(cl, *path, *sourcelist);
                  if (mergeResult < 0) {
                     Error("MergeRecursive", "error merging RNTuples");
                     return kFALSE;
                  }
                  return kTRUE;
               }
               TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
               Error("MergeRecursive", "Merging objects that don't inherit from TObject is unimplemented (key: %s of type %s in file %s)",
                      key->GetName(), key->GetClassName(), nextsource->GetName());
               canBeMerged = kFALSE;
            } else if (cl->IsTObject() && cl->GetMerge()) {
               // Check if already treated
               if (alreadyseen) continue;

               TList inputs;
               Bool_t oneGo = fHistoOneGo && cl->InheritsFrom(R__TH1_Class);

               // Loop over all source files and merge same-name object
               TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
               if (nextsource == 0) {
                  // There is only one file in the list
                  ROOT::MergeFunc_t func = cl->GetMerge();
                  func(obj, &inputs, &info);
                  info.fIsFirst = kFALSE;
               } else {
                  do {
                     // make sure we are at the correct directory level by cd'ing to path
                     TDirectory *ndir = nextsource->GetDirectory(path);
                     if (ndir) {
                        // For consistency (and persformance), we reset the MustCleanup be also for those
                        // 'key' retrieved indirectly.
                        // ndir->ResetBit(kMustCleanup);
                        ndir->cd();
                        TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
                        if (key2) {
                           TObject *hobj = key2->ReadObj();
                           if (!hobj) {
                              Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
                                   key->GetName(), key->GetTitle(), nextsource->GetName());
                              nextsource = (TFile*)sourcelist->After(nextsource);
                              continue;
                           }
                           // Set ownership for collections
                           if (hobj->InheritsFrom(TCollection::Class())) {
                              ((TCollection*)hobj)->SetOwner();
                           }
                           hobj->ResetBit(kMustCleanup);
                           inputs.Add(hobj);
                           if (!oneGo) {
                              ROOT::MergeFunc_t func = cl->GetMerge();
                              Long64_t result = func(obj, &inputs, &info);
                              info.fIsFirst = kFALSE;
                              if (result < 0) {
                                 Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
                                       key->GetName(), nextsource->GetName());
                              }
                              inputs.Delete();
                           }
                        }
                     }
                     nextsource = (TFile*)sourcelist->After( nextsource );
                  } while (nextsource);
                  // Merge the list, if still to be done
                  if (oneGo || info.fIsFirst) {
                     ROOT::MergeFunc_t func = cl->GetMerge();
                     func(obj, &inputs, &info);
                     info.fIsFirst = kFALSE;
                     inputs.Delete();
                  }
               }
            } else if (cl->IsTObject()) {
               // try synthesizing the Merge method call according to the TObject
               TList listH;
               TString listHargs;
               if (cl->GetMethodWithPrototype("Merge", "TCollection*,TFileMergeInfo*")) {
                  listHargs.Form("(TCollection*)0x%lx,(TFileMergeInfo*)0x%lx",
                                  (ULong_t)&listH, (ULong_t)&info);
               } else if (cl->GetMethodWithPrototype("Merge", "TCollection*")) {
                  listHargs.Form("((TCollection*)0x%lx)", (ULong_t)&listH);
               } else {
                  // pass unmergeable objects through to the output file
                  canBeMerged = kFALSE;
               }
               if (canBeMerged) {
                  if (alreadyseen) {
                     // skip already seen mergeable objects, don't skip unmergeable objects
                     continue;
                  }
                  // Loop over all source files and merge same-name object
                  TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
                  if (nextsource == 0) {
                     // There is only one file in the list
                     Int_t error = 0;
                     obj->Execute("Merge", listHargs.Data(), &error);
                     info.fIsFirst = kFALSE;
                     if (error) {
                        Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
                              obj->GetName(), key->GetName());
                     }
                  } else {
                     while (nextsource) {
                        // make sure we are at the correct directory level by cd'ing to path
                        TDirectory *ndir = nextsource->GetDirectory(path);
                        if (ndir) {
                           ndir->cd();
                           TKey *key2 = (TKey*)ndir->GetListOfKeys()->FindObject(key->GetName());
                           if (key2) {
                              TObject *hobj = key2->ReadObj();
                              if (!hobj) {
                                 Info("MergeRecursive", "could not read object for key {%s, %s}; skipping file %s",
                                      key->GetName(), key->GetTitle(), nextsource->GetName());
                                 nextsource = (TFile*)sourcelist->After(nextsource);
                                 continue;
                              }
                              // Set ownership for collections
                              if (hobj->InheritsFrom(TCollection::Class())) {
                                 ((TCollection*)hobj)->SetOwner();
                              }
                              hobj->ResetBit(kMustCleanup);
                              listH.Add(hobj);
                              Int_t error = 0;
                              obj->Execute("Merge", listHargs.Data(), &error);
                              info.fIsFirst = kFALSE;
                              if (error) {
                                 Error("MergeRecursive", "calling Merge() on '%s' with the corresponding object in '%s'",
                                       obj->GetName(), nextsource->GetName());
                              }
                              listH.Delete();
                           }
                        }
                        nextsource = (TFile*)sourcelist->After( nextsource );
                     }
                     // Merge the list, if still to be done
                     if (info.fIsFirst) {
                        Int_t error = 0;
                        obj->Execute("Merge", listHargs.Data(), &error);
                        info.fIsFirst = kFALSE;
                        listH.Delete();
                     }
                  }
               }
            } else {
               // Object is of no type that we can merge
               canBeMerged = kFALSE;
            }

            // now write the merged histogram (which is "in" obj) to the target file
            // note that this will just store obj in the current directory level,
            // which is not persistent until the complete directory itself is stored
            // by "target->SaveSelf()" below
            target->cd();

            oldkeyname = key->GetName();
            //!!if the object is a tree, it is stored in globChain...
            if(cl->InheritsFrom( TDirectory::Class() )) {
               //printf("cas d'une directory\n");

               auto dirobj = dynamic_cast<TDirectory*>(obj);
               TString dirpath(dirobj->GetPath());
               // coverity[unchecked_value] 'target' is from a file so GetPath always returns path starting with filename:
               dirpath.Remove(0, std::strlen(dirobj->GetFile()->GetPath()));

               // Do not delete the directory if it is part of the output
               // and we are in incremental mode (because it will be reuse
               // and has not been written to disk (for performance reason).
               // coverity[var_deref_model] the IsA()->InheritsFrom guarantees that the dynamic_cast will succeed.
               if (!(type&kIncremental) || dirobj->GetFile() != target) {
                  dirobj->ResetBit(kMustCleanup);
                  delete dirobj;
               }
               // Let's also delete the directory from the other source (thanks to the 'allNames'
               // mechanism above we will not process the directories when tranversing the next
               // files).
               TFile *nextsource = current_file ? (TFile*)sourcelist->After( current_file ) : (TFile*)sourcelist->First();
               while (nextsource) {
                  TDirectory *ndir = nextsource->GetDirectory(dirpath);
                  // For consistency (and persformance), we reset the MustCleanup be also for those
                  // 'key' retrieved indirectly.
                  if (ndir) {
                     ndir->ResetBit(kMustCleanup);
                     delete ndir;
                  }
                  nextsource = (TFile*)sourcelist->After( nextsource );
               }
            } else if (cl->InheritsFrom( TCollection::Class() )) {
               // Don't overwrite, if the object were not merged.
               if ( obj->Write( oldkeyname, canBeMerged ? TObject::kSingleKey | TObject::kOverwrite : TObject::kSingleKey) <= 0 ) {
                  status = kFALSE;
               }
               ((TCollection*)obj)->SetOwner();
               delete obj;
            } else {
               // Don't overwrite, if the object were not merged.
               // NOTE: this is probably wrong for emulated objects.
               if (cl->IsTObject()) {
                  if ( obj->Write( oldkeyname, canBeMerged ? TObject::kOverwrite : 0) <= 0) {
                     status = kFALSE;
                  }
                  obj->ResetBit(kMustCleanup);
               } else {
                  if ( target->WriteObjectAny( (void*)obj, cl, oldkeyname, canBeMerged ? "OverWrite" : "" ) <= 0) {
                     status = kFALSE;
                  }
               }
               cl->Destructor(obj); // just in case the class is not loaded.
            }
            info.Reset();
         } // while ( ( TKey *key = (TKey*)nextkey() ) )
      }
      current_file = current_file ? (TFile*)sourcelist->After(current_file) : (TFile*)sourcelist->First();
      if (current_file) {
         current_sourcedir = current_file->GetDirectory(path);
      } else {
         current_sourcedir = 0;
      }
   }
   // save modifications to the target directory.
   if (!(type&kIncremental)) {
      // In case of incremental build, we will call Write on the top directory/file, so we do not need
      // to call SaveSelf explicilty.
      target->SaveSelf(kTRUE);
   }

   return status;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge the files. If no output file was specified it will write into
/// the file "FileMerger.root" in the working directory. Returns true
/// on success, false in case of error.
/// The type is defined by the bit values in EPartialMergeType:
///   kRegular      : normal merge, overwritting the output file
///   kIncremental  : merge the input file with the content of the output file (if already exising) (default)
///   kAll          : merge all type of objects (default)
///   kResetable    : merge only the objects with a MergeAfterReset member function.
///   kNonResetable : merge only the objects without a MergeAfterReset member function.
///
/// If the type is set to kIncremental the output file is done deleted at the end of
/// this operation.  If the type is not set to kIncremental, the output file is closed.

Bool_t TFileMerger::PartialMerge(Int_t in_type)
{
   if (!fOutputFile) {
      TString outf(fOutputFilename);
      if (outf.IsNull()) {
         outf.Form("file:%s/FileMerger.root", gSystem->TempDirectory());
         Info("PartialMerge", "will merge the results to the file %s\n"
              "since you didn't specify a merge filename",
              TUrl(outf).GetFile());
      }
      if (!OutputFile(outf.Data())) {
         return kFALSE;
      }
   }

   // Special treament for the single file case ...
   if ((fFileList.GetEntries() == 1) && !fExcessFiles.GetEntries() &&
      !(in_type & kIncremental) && !fCompressionChange && !fExplicitCompLevel) {
      fOutputFile->Close();
      SafeDelete(fOutputFile);

      TFile *file = (TFile *) fFileList.First();
      if (!file || (file && file->IsZombie())) {
         Error("PartialMerge", "one-file case: problem attaching to file");
         return kFALSE;
      }
      Bool_t result = kTRUE;
      if (!(result = file->Cp(fOutputFilename))) {
         Error("PartialMerge", "one-file case: could not copy '%s' to '%s'",
                               file->GetPath(), fOutputFilename.Data());
         return kFALSE;
      }
      if (file->TestBit(kCanDelete)) file->Close();

      // Remove the temporary file
      if (fLocal && !file->InheritsFrom(TMemFile::Class())) {
         TUrl u(file->GetPath(), kTRUE);
         if (gSystem->Unlink(u.GetFile()) != 0)
            Warning("PartialMerge", "problems removing temporary local file '%s'", u.GetFile());
      }
      fFileList.Clear();
      return result;
   }

   fOutputFile->SetBit(kMustCleanup);

   TDirectory::TContext ctxt;

   Bool_t result = kTRUE;
   Int_t type = in_type;
   while (result && fFileList.GetEntries()>0) {
      result = MergeRecursive(fOutputFile, &fFileList, type);

      // Remove local copies if there are any
      TIter next(&fFileList);
      TFile *file;
      while ((file = (TFile*) next())) {
         // close the files
         if (file->TestBit(kCanDelete)) file->Close();
         // remove the temporary files
         if(fLocal && !file->InheritsFrom(TMemFile::Class())) {
            TString p(file->GetPath());
            // coverity[unchecked_value] Index is return a value with range or NPos to select the whole name.
            p = p(0, p.Index(':',0));
            gSystem->Unlink(p);
         }
      }
      fFileList.Clear();
      if (result && fExcessFiles.GetEntries() > 0) {
         // We merge the first set of files in the output,
         // we now need to open the next set and make
         // sure we accumulate into the output, so we
         // switch to incremental merging (if not already set)
         type = type | kIncremental;
         result = OpenExcessFiles();
      }
   }
   if (!result) {
      Error("Merge", "error during merge of your ROOT files");
   } else {
      // Close or write is required so the file is complete.
      if (in_type & kIncremental) {
         fOutputFile->Write("",TObject::kOverwrite);
      } else {
         gROOT->GetListOfFiles()->Remove(fOutputFile);
         fOutputFile->Close();
      }
   }

   // Cleanup
   if (in_type & kIncremental) {
      Clear();
   } else {
      fOutputFile->ResetBit(kMustCleanup);
      SafeDelete(fOutputFile);
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Open up to fMaxOpenedFiles of the excess files.

Bool_t TFileMerger::OpenExcessFiles()
{
   if (fPrintLevel > 0) {
      Printf("%s Opening the next %d files", fMsgPrefix.Data(), TMath::Min(fExcessFiles.GetEntries(), fMaxOpenedFiles - 1));
   }
   Int_t nfiles = 0;
   TIter next(&fExcessFiles);
   TObjString *url = 0;
   TString localcopy;
   // We want gDirectory untouched by anything going on here
   TDirectory::TContext ctxt;
   while( nfiles < (fMaxOpenedFiles-1) && ( url = (TObjString*)next() ) ) {
      TFile *newfile = 0;
      if (fLocal) {
         TUUID uuid;
         localcopy.Form("file:%s/ROOTMERGE-%s.root", gSystem->TempDirectory(), uuid.AsString());
         if (!TFile::Cp(url->GetName(), localcopy, url->TestBit(kCpProgress))) {
            Error("OpenExcessFiles", "cannot get a local copy of file %s", url->GetName());
            return kFALSE;
         }
         newfile = TFile::Open(localcopy, "READ");
      } else {
         newfile = TFile::Open(url->GetName(), "READ");
      }

      if (!newfile) {
         if (fLocal)
            Error("OpenExcessFiles", "cannot open local copy %s of URL %s",
                  localcopy.Data(), url->GetName());
         else
            Error("OpenExcessFiles", "cannot open file %s", url->GetName());
         return kFALSE;
      } else {
         if (fOutputFile && fOutputFile->GetCompressionLevel() != newfile->GetCompressionLevel()) fCompressionChange = kTRUE;

         newfile->SetBit(kCanDelete);
         fFileList.Add(newfile);
         ++nfiles;
         fExcessFiles.Remove(url);
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Intercept the case where the output TFile is deleted!

void TFileMerger::RecursiveRemove(TObject *obj)
{
   if (obj == fOutputFile) {
      Fatal("RecursiveRemove","Output file of the TFile Merger (targeting %s) has been deleted (likely due to a TTree larger than 100Gb)", fOutputFilename.Data());
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Set a limit to the number of files that TFileMerger will open simultaneously.
///
/// If the request is higher than the system limit, we reset it to the system limit.
/// If the request is less than two, we reset it to 2 (one for the output file and one for the input file).

void TFileMerger::SetMaxOpenedFiles(Int_t newmax)
{
   Int_t sysmax = R__GetSystemMaxOpenedFiles();
   if (newmax < sysmax) {
      fMaxOpenedFiles = newmax;
   } else {
      fMaxOpenedFiles = sysmax;
   }
   if (fMaxOpenedFiles < 2) {
      fMaxOpenedFiles = 2;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the prefix to be used when printing informational message.

void TFileMerger::SetMsgPrefix(const char *prefix)
{
   fMsgPrefix = prefix;
}

