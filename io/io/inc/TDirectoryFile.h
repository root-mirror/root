// @(#)root/io:$Id$
// Author: Rene Brun   22/01/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TDirectoryFile
#define ROOT_TDirectoryFile


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TDirectoryFile                                                       //
//                                                                      //
// Describe directory structure in a ROOT file.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Compression.h"
#include "TDirectory.h"
#include "TDatime.h"
#include "TList.h"

class TKey;
class TFile;

class TDirectoryFile : public TDirectory {

protected:
   Bool_t      fModified{kFALSE};        ///< True if directory has been modified
   Bool_t      fWritable{kFALSE};        ///< True if directory is writable
   TDatime     fDatimeC;                 ///< Date and time when directory is created
   TDatime     fDatimeM;                 ///< Date and time of last modification
   Int_t       fNbytesKeys{0};           ///< Number of bytes for the keys
   Int_t       fNbytesName{0};           ///< Number of bytes in TNamed at creation time
   Int_t       fBufferSize{0};           ///< Default buffer size to create new TKeys
   Long64_t    fSeekDir{0};              ///< Location of directory on file
   Long64_t    fSeekParent{0};           ///< Location of parent directory on file
   Long64_t    fSeekKeys{0};             ///< Location of Keys record on file
   TFile      *fFile{nullptr};           ///< Pointer to current file in memory
   TList      *fKeys{nullptr};           ///< Pointer to keys list in memory

   void        CleanTargets();
   void        InitDirectoryFile(TClass *cl = nullptr);
   void        BuildDirectoryFile(TFile* motherFile, TDirectory* motherDir);

private:
   TDirectoryFile(const TDirectoryFile &directory) = delete;  //Directories cannot be copied
   void operator=(const TDirectoryFile &) = delete; //Directories cannot be copied

public:
   // TDirectory status bits
   enum EStatusBits { kCloseDirectory = BIT(7) }; // Unused in ROOT, never set. Maybe only in external code.

   TDirectoryFile();
   TDirectoryFile(const char *name, const char *title, Option_t *option="", TDirectory* motherDir = nullptr);
   virtual ~TDirectoryFile();

          void        Append(TObject *obj, Bool_t replace = kFALSE) override;
          void        Add(TObject *obj, Bool_t replace = kFALSE) override { Append(obj,replace); }
          Int_t       AppendKey(TKey *key) override;
          void        Browse(TBrowser *b) override;
          void        Build(TFile* motherFile = nullptr, TDirectory* motherDir = nullptr) override { BuildDirectoryFile(motherFile, motherDir); }
          TObject    *CloneObject(const TObject *obj, Bool_t autoadd = kTRUE) override;
          void        Close(Option_t *option="") override;
          void        Copy(TObject &) const override { MayNotUse("Copy(TObject &)"); }
          Bool_t      cd(const char *path = nullptr) override;
          void        Delete(const char *namecycle="") override;
          void        FillBuffer(char *&buffer) override;
          TKey       *FindKey(const char *keyname) const override;
          TKey       *FindKeyAny(const char *keyname) const override;
          TObject    *FindObjectAny(const char *name) const override;
          TObject    *FindObjectAnyFile(const char *name) const override;
          TObject    *Get(const char *namecycle) override;
   /// See documentation of TDirectoryFile::Get(const char *namecycle)
   template <class T> inline T* Get(const char* namecycle)
   {
      return TDirectory::Get<T>(namecycle);
   }
           TDirectory *GetDirectory(const char *apath, Bool_t printError = false, const char *funcname = "GetDirectory") override;
           void       *GetObjectChecked(const char *namecycle, const char* classname) override;
           void       *GetObjectChecked(const char *namecycle, const TClass* cl) override;
           void       *GetObjectUnchecked(const char *namecycle) override;
           Int_t       GetBufferSize() const override;
   const TDatime      &GetCreationDate() const { return fDatimeC; }
           TFile      *GetFile() const override { return fFile; }
           TKey       *GetKey(const char *name, Short_t cycle=9999) const override;
           TList      *GetListOfKeys() const override { return fKeys; }
   const TDatime      &GetModificationDate() const { return fDatimeM; }
           Int_t       GetNbytesKeys() const override { return fNbytesKeys; }
           Int_t       GetNkeys() const override { return fKeys->GetSize(); }
           Long64_t    GetSeekDir() const override { return fSeekDir; }
           Long64_t    GetSeekParent() const override { return fSeekParent; }
           Long64_t    GetSeekKeys() const override { return fSeekKeys; }
           Bool_t      IsModified() const override { return fModified; }
           Bool_t      IsWritable() const override { return fWritable; }
           void        ls(Option_t *option="") const override;
           TDirectory *mkdir(const char *name, const char *title="", Bool_t returnExistingDirectory = kFALSE) override;
           TFile      *OpenFile(const char *name, Option_t *option= "",
                            const char *ftitle = "", Int_t compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault,
                            Int_t netopt = 0) override;
           void        Purge(Short_t nkeep=1) override;
           void        ReadAll(Option_t *option="") override;
           Int_t       ReadKeys(Bool_t forceRead=kTRUE) override;
           Int_t       ReadTObject(TObject *obj, const char *keyname) override;
   virtual void        ResetAfterMerge(TFileMergeInfo *);
           void        rmdir(const char *name) override;
           void        Save() override;
           void        SaveSelf(Bool_t force = kFALSE) override;
           Int_t       SaveObjectAs(const TObject *obj, const char *filename="", Option_t *option="") const override;
           void        SetBufferSize(Int_t bufsize) override;
           void        SetModified() override {fModified = kTRUE;}
           void        SetSeekDir(Long64_t v) override { fSeekDir = v; }
           void        SetTRefAction(TObject *ref, TObject *parent) override;
           void        SetWritable(Bool_t writable=kTRUE) override;
           Int_t       Sizeof() const override;
           Int_t       Write(const char *name=nullptr, Int_t opt=0, Int_t bufsize=0) override;
           Int_t       Write(const char *name=nullptr, Int_t opt=0, Int_t bufsize=0) const override;
           Int_t       WriteTObject(const TObject *obj, const char *name=nullptr, Option_t *option="", Int_t bufsize=0) override;
           Int_t       WriteObjectAny(const void *obj, const char *classname, const char *name, Option_t *option="", Int_t bufsize=0) override;
           Int_t       WriteObjectAny(const void *obj, const TClass *cl, const char *name, Option_t *option="", Int_t bufsize=0) override;
           void        WriteDirHeader() override;
           void        WriteKeys() override;

   ClassDefOverride(TDirectoryFile,5)  //Describe directory structure in a ROOT file
};

#endif
