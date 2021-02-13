// @(#)root/base:$Id$
// Author: Andreas-Joachim Peters   20/9/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFileInfo
#define ROOT_TFileInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFileInfo                                                            //
//                                                                      //
// Class describing a generic file including meta information.          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TList.h"

#ifdef R__LESS_INCLUDES
class TUrl;
class TUUID;
class TMD5;
#else
#include "TUrl.h"
#include "TUUID.h"
#include "TMD5.h"
#endif

class TFileInfoMeta;


class TFileInfo : public TNamed {

private:
   TUrl            *fCurrentUrl;     //! current URL to access the file, points to URL
                                     //  in the fUrlList or 0, if the list end is reached
   TList           *fUrlList;        // list of file URLs
   Long64_t         fSize;           // file size
   TUUID           *fUUID;           //-> uuid of the referenced file
   TMD5            *fMD5;            //-> md5 digest of the file
   TList           *fMetaDataList;   // generic list of file meta data object(s)

   Int_t            fIndex;          // Index to be used when sorting with index

   void             ParseInput(const char *in);

   TFileInfo& operator=(const TFileInfo&) = delete;

public:
   enum EStatusBits {
      kStaged    = BIT(15),
      kCorrupted = BIT(16),
      kSortWithIndex  = BIT(17)     // Use index when sorting (in Compare)
   };

   TFileInfo(const char *url = 0, Long64_t size = -1, const char *uuid = 0,
             const char *md5 = 0, TObject *meta = 0);
   TFileInfo(const TFileInfo &);

   virtual ~TFileInfo();

   void            ResetUrl() { fCurrentUrl = (TUrl*)fUrlList->First(); }
   TUrl           *NextUrl();
   TUrl           *GetCurrentUrl() const;
   TUrl           *GetFirstUrl() const { return (TUrl*)fUrlList->First(); }
   TUrl           *GetUrlAt(Int_t i) const { return (TUrl*)fUrlList->At(i); }
   Bool_t          RemoveUrlAt(Int_t i);
   Int_t           GetNUrls() const    { return fUrlList->GetEntries(); }

   Bool_t          SetCurrentUrl(const char *url);
   Bool_t          SetCurrentUrl(TUrl *url);

   Long64_t        GetSize() const         { return fSize; }
   TUUID          *GetUUID() const         { return fUUID; }
   TMD5           *GetMD5() const          { return fMD5; }
   TList          *GetMetaDataList() const { return fMetaDataList; }
   TFileInfoMeta  *GetMetaData(const char *meta = 0) const;

   void            SetSize(Long64_t size)  { fSize = size; }
   void            SetUUID(const char *uuid);

   TUrl           *FindByUrl(const char *url, Bool_t withDeflt = kFALSE);

   Bool_t          AddUrl(const char *url, Bool_t infront = kFALSE);
   Bool_t          RemoveUrl(const char *url);
   Bool_t          AddMetaData(TObject *meta);
   Bool_t          RemoveMetaData(const char *meta = 0);

   Bool_t          IsSortable() const { return kTRUE; }
   Int_t           Compare(const TObject *obj) const;

   Int_t           GetIndex() const { return fIndex; }
   void            SetIndex(Int_t idx) { fIndex = idx; }

   void            Print(Option_t *options="") const;

   ClassDef(TFileInfo,4)   // Describes generic file info including meta data information
};


class TFileInfoMeta : public TNamed {

private:
   Long64_t      fEntries;    // number of entries in tree or number of objects
   Long64_t      fFirst;      // first valid tree entry
   Long64_t      fLast;       // last valid tree entry
   Bool_t        fIsTree;     // true if type is a TTree (or TTree derived)
   Long64_t      fTotBytes;   // uncompressed size in bytes
   Long64_t      fZipBytes;   // compressed size in bytes

   TFileInfoMeta& operator=(const TFileInfoMeta &) = delete;

public:
   enum EStatusBits { kExternal  = BIT(15) };

   TFileInfoMeta() : fEntries(-1), fFirst(0), fLast(-1),
                     fIsTree(kFALSE), fTotBytes(-1), fZipBytes(-1)
                     { ResetBit(TFileInfoMeta::kExternal); }
   TFileInfoMeta(const char *objPath, const char *objClass = "TTree",
                 Long64_t entries = -1, Long64_t first = 0, Long64_t last = -1,
                 Long64_t totbytes = -1, Long64_t zipbytes = -1);
   TFileInfoMeta(const char *objPath, const char *objDir,
                 const char *objClass, Long64_t entries = -1,
                 Long64_t first = 0, Long64_t last = -1,
                 Long64_t totbytes = -1, Long64_t zipbytes = -1);
   TFileInfoMeta(const TFileInfoMeta &m);

   virtual ~TFileInfoMeta() { }

   const char     *GetObject() const;
   const char     *GetClass() const        { return GetTitle(); }
   const char     *GetDirectory() const;
   Long64_t        GetEntries() const      { return fEntries; }
   Long64_t        GetFirst() const        { return fFirst; }
   Long64_t        GetLast() const         { return fLast; }
   Bool_t          IsTree() const          { return fIsTree; }
   Long64_t        GetTotBytes() const     { return fTotBytes; }
   Long64_t        GetZipBytes() const     { return fZipBytes; }

   void            SetEntries(Long64_t entries) { fEntries = entries; }
   void            SetFirst(Long64_t first)     { fFirst = first; }
   void            SetLast(Long64_t last)       { fLast = last; }
   void            SetTotBytes(Long64_t tot)    { fTotBytes = tot; }
   void            SetZipBytes(Long64_t zip)    { fZipBytes = zip; }

   void            Print(Option_t *options="") const;

   ClassDef(TFileInfoMeta,2)   // Describes TFileInfo meta data
};

#endif
