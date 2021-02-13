// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranch
#define ROOT_TBranch


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBranch                                                              //
//                                                                      //
// A TTree object is a list of TBranchs.                                //
//   A TBranch describes the branch data structure and supports :       //
//     the list of TBaskets (branch buffers) associated to this branch. //
//     the list of TLeaves (branch description)                         //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TAttFill.h"
#include "TObjArray.h"
#include "TBranchCacheInfo.h"
#include "TDataType.h"
#include "Compression.h"
#include "ROOT/TIOFeatures.hxx"

class TTree;
class TBasket;
class TBranchElement;
class TLeaf;
class TBrowser;
class TDirectory;
class TFile;
class TClonesArray;
class TTreeCloner;
class TTreeCache;

namespace ROOT {
namespace Experimental {
namespace Internal {
class TBulkBranchRead;
}
}
namespace Internal {
class TBranchIMTHelper; ///< A helper class for managing IMT work during TTree:Fill operations.
}
}

const Int_t kDoNotProcess = BIT(10); // Active bit for branches
const Int_t kIsClone      = BIT(11); // to indicate a TBranchClones
const Int_t kBranchObject = BIT(12); // branch is a TObject*
const Int_t kBranchAny    = BIT(17); // branch is an object*
const Int_t kMapObject    = kBranchObject | kBranchAny;

namespace ROOT {
namespace Experimental {
namespace Internal {

///\class TBulkBranchRead
/// Helper class for reading many branch entries at once to optimize throughput.
class TBulkBranchRead {

   friend class ::TBranch;

public:
   Int_t GetBulkEntries(Long64_t evt, TBuffer &user_buf);
   Int_t GetEntriesSerialized(Long64_t evt, TBuffer &user_buf);
   Int_t GetEntriesSerialized(Long64_t evt, TBuffer &user_buf, TBuffer *count_buf);
   Bool_t SupportsBulkRead() const;

private:
   TBulkBranchRead(TBranch &parent)
      : fParent(parent)
   {}

   TBranch &fParent;
};
}
}
}

class TBranch : public TNamed , public TAttFill {
   using TIOFeatures = ROOT::TIOFeatures;

protected:
   friend class TTreeCache;
   friend class TTreeCloner;
   friend class TTree;
   friend class TBranchElement;
   friend class ROOT::Experimental::Internal::TBulkBranchRead;

   // TBranch status bits
   enum EStatusBits {
      kDoNotProcess = ::kDoNotProcess, // Active bit for branches
      kIsClone      = ::kIsClone,      // to indicate a TBranchClones
      kBranchObject = ::kBranchObject, // branch is a TObject*
      kBranchAny    = ::kBranchAny,    // branch is an object*
      // kMapObject    = kBranchObject | kBranchAny;
      kAutoDelete   = BIT(15),

      kDoNotUseBufferMap = BIT(22) // If set, at least one of the entry in the branch will use the buffer's map of classname and objects.
   };

   using BulkObj = ROOT::Experimental::Internal::TBulkBranchRead;
   static Int_t fgCount;          ///<! branch counter
   Int_t       fCompress;         ///<  Compression level and algorithm
   Int_t       fBasketSize;       ///<  Initial Size of  Basket Buffer
   Int_t       fEntryOffsetLen;   ///<  Initial Length of fEntryOffset table in the basket buffers
   Int_t       fWriteBasket;      ///<  Last basket number written
   Long64_t    fEntryNumber;      ///<  Current entry number (last one filled in this branch)
   TBasket    *fExtraBasket;      ///<! Allocated basket not currently holding any data.
   TIOFeatures fIOFeatures;       ///<  IO features for newly-created baskets.
   Int_t       fOffset;           ///<  Offset of this branch
   Int_t       fMaxBaskets;       ///<  Maximum number of Baskets so far
   Int_t       fNBaskets;         ///<! Number of baskets in memory
   Int_t       fSplitLevel;       ///<  Branch split level
   Int_t       fNleaves;          ///<! Number of leaves
   Int_t       fReadBasket;       ///<! Current basket number when reading
   Long64_t    fReadEntry;        ///<! Current entry number when reading
   Long64_t    fFirstBasketEntry; ///<! First entry in the current basket.
   Long64_t    fNextBasketEntry;  ///<! Next entry that will requires us to go to the next basket
   TBasket    *fCurrentBasket;    ///<! Pointer to the current basket.
   Long64_t    fEntries;          ///<  Number of entries
   Long64_t    fFirstEntry;       ///<  Number of the first entry in this branch
   Long64_t    fTotBytes;         ///<  Total number of bytes in all leaves before compression
   Long64_t    fZipBytes;         ///<  Total number of bytes in all leaves after compression
   TObjArray   fBranches;         ///< -> List of Branches of this branch
   TObjArray   fLeaves;           ///< -> List of leaves of this branch
   TObjArray   fBaskets;          ///< -> List of baskets of this branch
   Int_t      *fBasketBytes;      ///<[fMaxBaskets] Length of baskets on file
   Long64_t   *fBasketEntry;      ///<[fMaxBaskets] Table of first entry in each basket
   Long64_t   *fBasketSeek;       ///<[fMaxBaskets] Addresses of baskets on file
   TTree      *fTree;             ///<! Pointer to Tree header
   TBranch    *fMother;           ///<! Pointer to top-level parent branch in the tree.
   TBranch    *fParent;           ///<! Pointer to parent branch.
   char       *fAddress;          ///<! Address of 1st leaf (variable or object)
   TDirectory *fDirectory;        ///<! Pointer to directory where this branch buffers are stored
   TString     fFileName;         ///<  Name of file where buffers are stored ("" if in same file as Tree header)
   TBuffer    *fEntryBuffer;      ///<! Buffer used to directly pass the content without streaming
   TBuffer    *fTransientBuffer;  ///<! Pointer to the current transient buffer.
   TList      *fBrowsables;       ///<! List of TVirtualBranchBrowsables used for Browse()
   BulkObj     fBulk;             ///<! Helper for performing bulk IO

   Bool_t      fSkipZip;          ///<! After being read, the buffer will not be unzipped.

   using CacheInfo_t = ROOT::Internal::TBranchCacheInfo;
   CacheInfo_t fCacheInfo;        ///<! Hold info about which basket are in the cache and if they have been retrieved from the cache.

   typedef void (TBranch::*ReadLeaves_t)(TBuffer &b);
   ReadLeaves_t fReadLeaves;      ///<! Pointer to the ReadLeaves implementation to use.
   typedef void (TBranch::*FillLeaves_t)(TBuffer &b);
   FillLeaves_t fFillLeaves;      ///<! Pointer to the FillLeaves implementation to use.
   void     ReadLeavesImpl(TBuffer &b);
   void     ReadLeaves0Impl(TBuffer &b);
   void     ReadLeaves1Impl(TBuffer &b);
   void     ReadLeaves2Impl(TBuffer &b);
   void     FillLeavesImpl(TBuffer &b);

   void     SetSkipZip(Bool_t skip = kTRUE) { fSkipZip = skip; }
   void     Init(const char *name, const char *leaflist, Int_t compress);

   TBasket *GetFreshBasket(Int_t basketnumber, TBuffer *user_buffer);
   TBasket *GetFreshCluster();
   Int_t    WriteBasket(TBasket* basket, Int_t where) { return WriteBasketImpl(basket, where, nullptr); }

   TString  GetRealFileName() const;

   virtual void SetAddressImpl(void *addr, Bool_t /* implied */) { SetAddress(addr); }

private:
   Int_t    GetBasketAndFirst(TBasket*& basket, Long64_t& first, TBuffer* user_buffer);
   TBasket *GetBasketImpl(Int_t basket, TBuffer* user_buffer);
   Int_t    GetBulkEntries(Long64_t, TBuffer&);
   Int_t    GetEntriesSerialized(Long64_t N, TBuffer& user_buf) {return GetEntriesSerialized(N, user_buf, nullptr);}
   Int_t    GetEntriesSerialized(Long64_t, TBuffer&, TBuffer*);
   Int_t    FillEntryBuffer(TBasket* basket,TBuffer* buf, Int_t& lnew);
   Int_t    WriteBasketImpl(TBasket* basket, Int_t where, ROOT::Internal::TBranchIMTHelper *);
   TBranch(const TBranch&) = delete;             // not implemented
   TBranch& operator=(const TBranch&) = delete;  // not implemented

public:
   TBranch();
   TBranch(TTree *tree, const char *name, void *address, const char *leaflist, Int_t basketsize=32000, Int_t compress = ROOT::RCompressionSetting::EAlgorithm::kInherit);
   TBranch(TBranch *parent, const char *name, void *address, const char *leaflist, Int_t basketsize=32000, Int_t compress = ROOT::RCompressionSetting::EAlgorithm::kInherit);
   virtual ~TBranch();

   virtual void      AddBasket(TBasket &b, Bool_t ondisk, Long64_t startEntry);
   virtual void      AddLastBasket(Long64_t startEntry);
           Int_t     BackFill();
   virtual void      Browse(TBrowser *b);
   virtual void      DeleteBaskets(Option_t* option="");
   virtual void      DropBaskets(Option_t *option = "");
           void      ExpandBasketArrays();
           Int_t     Fill() { return FillImpl(nullptr); }
   virtual Int_t     FillImpl(ROOT::Internal::TBranchIMTHelper *);
   virtual TBranch  *FindBranch(const char *name);
   virtual TLeaf    *FindLeaf(const char *name);
           Int_t     FlushBaskets();
           Int_t     FlushOneBasket(UInt_t which);

   virtual char     *GetAddress() const {return fAddress;}
           TBasket  *GetBasket(Int_t basket) {return GetBasketImpl(basket, nullptr);}
           Int_t    *GetBasketBytes() const {return fBasketBytes;}
           Long64_t *GetBasketEntry() const {return fBasketEntry;}
   virtual Long64_t  GetBasketSeek(Int_t basket) const;
   virtual Int_t     GetBasketSize() const {return fBasketSize;}
           ROOT::Experimental::Internal::TBulkBranchRead &GetBulkRead() { return fBulk; }
   virtual TList    *GetBrowsables();
   virtual const char* GetClassName() const;
           Int_t     GetCompressionAlgorithm() const;
           Int_t     GetCompressionLevel() const;
           Int_t     GetCompressionSettings() const;
   TDirectory       *GetDirectory() const {return fDirectory;}
   virtual Int_t     GetEntry(Long64_t entry=0, Int_t getall = 0);
   virtual Int_t     GetEntryExport(Long64_t entry, Int_t getall, TClonesArray *list, Int_t n);
           Int_t     GetEntryOffsetLen() const { return fEntryOffsetLen; }
           Int_t     GetEvent(Long64_t entry=0) {return GetEntry(entry);}
   virtual TString   GetFullName() const;
   const char       *GetIconName() const;
   virtual Int_t     GetExpectedType(TClass *&clptr,EDataType &type);
   virtual TLeaf    *GetLeaf(const char *name) const;
   virtual TFile    *GetFile(Int_t mode=0);
   const char       *GetFileName()    const {return fFileName.Data();}
           Int_t     GetOffset()      const {return fOffset;}
           Int_t     GetReadBasket()  const {return fReadBasket;}
           Long64_t  GetReadEntry()   const {return fReadEntry;}
           Int_t     GetWriteBasket() const {return fWriteBasket;}
           Long64_t  GetTotalSize(Option_t *option="")   const;
           Long64_t  GetTotBytes(Option_t *option="")    const;
           Long64_t  GetZipBytes(Option_t *option="")    const;
           Long64_t  GetEntryNumber() const {return fEntryNumber;}
           Long64_t  GetFirstEntry()  const {return fFirstEntry; }
         TIOFeatures GetIOFeatures() const;
         TObjArray  *GetListOfBaskets()  {return &fBaskets;}
         TObjArray  *GetListOfBranches() {return &fBranches;}
         TObjArray  *GetListOfLeaves()   {return &fLeaves;}
           Int_t     GetMaxBaskets()  const  {return fMaxBaskets;}
           Int_t     GetNleaves()     const {return fNleaves;}
           Int_t     GetSplitLevel()  const {return fSplitLevel;}
           Long64_t  GetEntries()     const {return fEntries;}
           TTree    *GetTree()        const {return fTree;}
   virtual Int_t     GetRow(Int_t row);
   virtual Bool_t    GetMakeClass() const;
   TBranch          *GetMother() const;
   TBranch          *GetSubBranch(const TBranch *br) const;
   TBuffer          *GetTransientBuffer(Int_t size);
   Bool_t            IsAutoDelete() const;
   Bool_t            IsFolder() const;
   virtual void      KeepCircular(Long64_t maxEntries);
   virtual Int_t     LoadBaskets();
   virtual void      Print(Option_t *option="") const;
           void      PrintCacheInfo() const;
   virtual void      ReadBasket(TBuffer &b);
   virtual void      Refresh(TBranch *b);
   virtual void      Reset(Option_t *option="");
   virtual void      ResetAfterMerge(TFileMergeInfo *);
   virtual void      ResetAddress();
   virtual void      ResetReadEntry() {fReadEntry = -1;}
   virtual void      SetAddress(void *add);
   virtual void      SetObject(void *objadd);
   virtual void      SetAutoDelete(Bool_t autodel=kTRUE);
   virtual void      SetBasketSize(Int_t buffsize);
   virtual void      SetBufferAddress(TBuffer *entryBuffer);
   void              SetCompressionAlgorithm(Int_t algorithm = ROOT::RCompressionSetting::EAlgorithm::kUseGlobal);
   void              SetCompressionLevel(Int_t level = ROOT::RCompressionSetting::ELevel::kUseMin);
   void              SetCompressionSettings(Int_t settings = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault);
   virtual void      SetEntries(Long64_t entries);
   virtual void      SetEntryOffsetLen(Int_t len, Bool_t updateSubBranches = kFALSE);
   virtual void      SetFirstEntry( Long64_t entry );
   virtual void      SetFile(TFile *file=0);
   virtual void      SetFile(const char *filename);
   void              SetIOFeatures(TIOFeatures &features) {fIOFeatures = features;}
   virtual Bool_t    SetMakeClass(Bool_t decomposeObj = kTRUE);
   virtual void      SetOffset(Int_t offset=0) {fOffset=offset;}
   virtual void      SetStatus(Bool_t status=1);
   virtual void      SetTree(TTree *tree) { fTree = tree;}
   virtual void      SetupAddresses();
           Bool_t    SupportsBulkRead() const;
   virtual void      UpdateAddress() {;}
   virtual void      UpdateFile();

   static  void      ResetCount();

   ClassDef(TBranch, 13); // Branch descriptor
};

//______________________________________________________________________________
inline Int_t TBranch::GetCompressionAlgorithm() const
{
   return (fCompress < 0) ? -1 : fCompress / 100;
}

//______________________________________________________________________________
inline Int_t TBranch::GetCompressionLevel() const
{
   return (fCompress < 0) ? -1 : fCompress % 100;
}

//______________________________________________________________________________
inline Int_t TBranch::GetCompressionSettings() const
{
   return (fCompress < 0) ? -1 : fCompress;
}

namespace ROOT {
namespace Experimental {
namespace Internal {

inline Int_t  TBulkBranchRead::GetBulkEntries(Long64_t evt, TBuffer& user_buf) { return fParent.GetBulkEntries(evt, user_buf); }
inline Int_t  TBulkBranchRead::GetEntriesSerialized(Long64_t evt, TBuffer& user_buf) { return fParent.GetEntriesSerialized(evt, user_buf); }
inline Int_t  TBulkBranchRead::GetEntriesSerialized(Long64_t evt, TBuffer& user_buf, TBuffer* count_buf) { return fParent.GetEntriesSerialized(evt, user_buf, count_buf); }
inline Bool_t TBulkBranchRead::SupportsBulkRead() const { return fParent.SupportsBulkRead(); }

}  // Internal
}  // Experimental
}  // ROOT

#endif
