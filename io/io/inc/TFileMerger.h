// @(#)root/io:$Id$
// Author: Andreas Peters + Fons Rademakers   26/5/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFileMerger
#define ROOT_TFileMerger

#include "TObject.h"
#include "TList.h"
#include "TString.h"
#include "TStopwatch.h"
#include <string>

class TFile;
class TDirectory;

namespace ROOT {
class TIOFeatures;
}  // namespace ROOT

class TFileMerger : public TObject {
private:
   using TIOFeatures = ROOT::TIOFeatures;

   TFileMerger(const TFileMerger&) = delete;
   TFileMerger& operator=(const TFileMerger&) = delete;

protected:
   TStopwatch     fWatch;                     ///< Stop watch to measure file copy speed
   TList          fFileList;                  ///< A list the file (TFile*) which shall be merged
   TFile         *fOutputFile{nullptr};       ///< The outputfile for merging
   TString        fOutputFilename;            ///< The name of the outputfile for merging
   Bool_t         fFastMethod{kTRUE};         ///< True if using Fast merging algorithm (default)
   Bool_t         fNoTrees{kFALSE};           ///< True if Trees should not be merged (default is kFALSE)
   Bool_t         fExplicitCompLevel{kFALSE}; ///< True if the user explicitly requested a compressio level change (default kFALSE)
   Bool_t         fCompressionChange{kFALSE}; ///< True if the output and input have different compression level (default kFALSE)
   Int_t          fPrintLevel{0};             ///< How much information to print out at run time
   TString        fMergeOptions;              ///< Options (in string format) to be passed down to the Merge functions
   TIOFeatures   *fIOFeatures{nullptr};       ///< IO features to use in the output file.
   TString        fMsgPrefix{"TFileMerger"};  ///< Prefix to be used when printing informational message (default TFileMerger)

   Int_t          fMaxOpenedFiles;            ///< Maximum number of files opened at the same time by the TFileMerger
   Bool_t         fLocal;                     ///< Makes local copies of merging files if True (default is kTRUE)
   Bool_t         fHistoOneGo;                ///< Merger histos in one go (default is kTRUE)
   TString        fObjectNames;               ///< List of object names to be either merged exclusively or skipped
   TList          fMergeList;                 ///< list of TObjString containing the name of the files need to be merged
   TList          fExcessFiles;               ///<! List of TObjString containing the name of the files not yet added to fFileList due to user or system limitiation on the max number of files opened.

   Bool_t         OpenExcessFiles();
   virtual Bool_t AddFile(TFile *source, Bool_t own, Bool_t cpProgress);
   virtual Bool_t MergeRecursive(TDirectory *target, TList *sourcelist, Int_t type = kRegular | kAll);

public:
   /// Type of the partial merge
   enum EPartialMergeType {
      kRegular      = 0,             ///< Normal merge, overwritting the output file.
      kIncremental  = BIT(1),        ///< Merge the input file with the content of the output file (if already exising).
      kResetable    = BIT(2),        ///< Only the objects with a MergeAfterReset member function.
      kNonResetable = BIT(3),        ///< Only the objects without a MergeAfterReset member function.

      kAll            = BIT(2)|BIT(3),      ///< Merge all type of objects (default)
      kAllIncremental = kIncremental | kAll, ///< Merge incrementally all type of objects.

      kOnlyListed     = BIT(4),        ///< Only the objects specified in fObjectNames list
      kSkipListed     = BIT(5),        ///< Skip objects specified in fObjectNames list
      kKeepCompression= BIT(6)         ///< Keep compression level unchanged for each input files
   };

   TFileMerger(Bool_t isLocal = kTRUE, Bool_t histoOneGo = kTRUE);
   virtual ~TFileMerger();

   Int_t       GetPrintLevel() const { return fPrintLevel; }
   void        SetPrintLevel(Int_t level) { fPrintLevel = level; }
   Bool_t      HasCompressionChange() const { return fCompressionChange; }
   const char *GetOutputFileName() const { return fOutputFilename; }
   TList      *GetMergeList() { return &fMergeList; }
   TFile      *GetOutputFile() const { return fOutputFile; }
   Int_t       GetMaxOpenedFiles() const { return fMaxOpenedFiles; }
   void        SetMaxOpenedFiles(Int_t newmax);
   const char *GetMsgPrefix() const { return fMsgPrefix; }
   void        SetMsgPrefix(const char *prefix);
   const char *GetMergeOptions() { return fMergeOptions; }
   void        SetMergeOptions(const TString &options) { fMergeOptions = options; }
   void        SetMergeOptions(const std::string_view &options) { fMergeOptions = options; }
   void        SetIOFeatures(ROOT::TIOFeatures &features) { fIOFeatures = &features; }
   void        AddObjectNames(const char *name) {fObjectNames += name; fObjectNames += " ";}
   const char *GetObjectNames() const {return fObjectNames.Data();}
   void        ClearObjectNames() {fObjectNames.Clear();}

    //--- file management interface
   virtual Bool_t SetCWD(const char * /*path*/) { MayNotUse("SetCWD"); return kFALSE; }
   virtual const char *GetCWD() { MayNotUse("GetCWD"); return 0; }

   //--- file merging interface
   virtual void   Reset();
   virtual Bool_t AddFile(const char *url, Bool_t cpProgress = kTRUE);
   virtual Bool_t AddFile(TFile *source, Bool_t cpProgress = kTRUE);
   virtual Bool_t AddAdoptFile(TFile *source, Bool_t cpProgress = kTRUE);
   virtual Bool_t OutputFile(const char *url, Bool_t force);
   virtual Bool_t OutputFile(const char *url, Bool_t force, Int_t compressionLevel);
   virtual Bool_t OutputFile(const char *url, const char *mode = "RECREATE");
   virtual Bool_t OutputFile(const char *url, const char *mode, Int_t compressionLevel);
   virtual Bool_t OutputFile(std::unique_ptr<TFile> file);
   virtual void   PrintFiles(Option_t *options);
   virtual Bool_t Merge(Bool_t = kTRUE);
   virtual Bool_t PartialMerge(Int_t type = kAll | kIncremental);
   virtual void   SetFastMethod(Bool_t fast=kTRUE)  {fFastMethod = fast;}
   virtual void   SetNotrees(Bool_t notrees=kFALSE) {fNoTrees = notrees;}
   virtual void        RecursiveRemove(TObject *obj);

   ClassDef(TFileMerger, 6)  // File copying and merging services
};

#endif

