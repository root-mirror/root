// @(#)root/base:$Id$
// Author: Jan Fiete Grosse-Oetringhaus  01/06/07

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFileCollection
#define ROOT_TFileCollection

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFileCollection                                                      //
//                                                                      //
// Class that contains a list of TFileInfo's and accumulated meta       //
// data information about its entries. This class is used to describe   //
// file sets as stored by Grid file catalogs, by PROOF or any other     //
// collection of TFile names.                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#include "TString.h"

class THashList;
class TMap;
class TList;
class TCollection;
class TFileInfo;
class TFileInfoMeta;
class TObjString;


class TFileCollection : public TNamed {

private:
   THashList  *fList;               //-> list of TFileInfos
   TList      *fMetaDataList;       //-> generic list of file meta data object(s)
                                    //  (summed over entries of fList)
   TString     fDefaultTree;        // name of default tree
   Long64_t    fTotalSize;          // total size of files in the list
   Long64_t    fNFiles;             // number of files ( == fList->GetEntries(), needed
                                    // because TFileCollection might be read without fList)
   Long64_t    fNStagedFiles;       // number of staged files
   Long64_t    fNCorruptFiles;      // number of corrupt files

   TFileCollection(const TFileCollection&) = delete;
   TFileCollection& operator=(const TFileCollection&) = delete;

   void PrintDetailed(TString &showOnly) const;
   void FormatSize(Long64_t bytes, TString &um, Double_t &size) const;

public:
   enum EStatusBits {
      kRemoteCollection = BIT(15)   // the collection is not staged
   };
   TFileCollection(const char *name = nullptr, const char *title = nullptr,
                   const char *file = nullptr, Int_t nfiles = -1, Int_t firstfile = 1);
   virtual ~TFileCollection();

   Int_t           Add(TFileInfo *info);
   Int_t           Add(TFileCollection *coll);
   Int_t           AddFromFile(const char *file, Int_t nfiles = -1, Int_t firstfile = 1);
   Int_t           Add(const char *path);
   THashList      *GetList() { return fList; }
   void            SetList(THashList* list) { fList = list; }

   TObjString     *ExportInfo(const char *name = 0, Int_t popt = 0);

   Long64_t        Merge(TCollection* list);
   Int_t           RemoveDuplicates();
   Int_t           Update(Long64_t avgsize = -1);
   void            Sort(Bool_t useindex = kFALSE);
   void            SetAnchor(const char *anchor);
   void            Print(Option_t *option = "") const;

   void            SetBitAll(UInt_t f);
   void            ResetBitAll(UInt_t f);

   Long64_t        GetTotalSize() const           { return fTotalSize; }
   Long64_t        GetNFiles() const              { return fNFiles; }
   Long64_t        GetNStagedFiles() const        { return fNStagedFiles; }
   Long64_t        GetNCorruptFiles() const       { return fNCorruptFiles; }
   Float_t         GetStagedPercentage() const
                   { return (fNFiles > 0) ? 100. * fNStagedFiles / fNFiles : 0; }
   Float_t         GetCorruptedPercentage() const
                   { return (fNFiles > 0) ? 100. * fNCorruptFiles / fNFiles : 0; }

   const char     *GetDefaultTreeName() const;
   void            SetDefaultTreeName(const char* treeName) { fDefaultTree = treeName; }
   Long64_t        GetTotalEntries(const char *tree) const;

   TFileInfoMeta  *GetMetaData(const char *meta = 0) const;
   void            SetDefaultMetaData(const char *meta);
   Bool_t          AddMetaData(TObject *meta);
   void            RemoveMetaData(const char *meta = 0);

   TFileCollection *GetStagedSubset();

   TFileCollection *GetFilesOnServer(const char *server);
   TMap            *GetFilesPerServer(const char *exclude = 0, Bool_t curronly =  kFALSE);

   ClassDef(TFileCollection, 3)  // Collection of TFileInfo objects
};

#endif
