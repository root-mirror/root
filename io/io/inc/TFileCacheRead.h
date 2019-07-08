/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFileCacheRead
#define ROOT_TFileCacheRead

#include "TObject.h"

#include "TFile.h"

#include <memory>

class TBranch;
class TFilePrefetch;
class TFileBufferRead;

class TFileCacheRead : public TObject {

protected:
   TFilePrefetch *fPrefetch;         ///<! Object that does the asynchronous reading in another thread
   Int_t          fBufferSizeMin;    ///< Original size of fBuffer
   Int_t          fBufferSize;       ///< Allocated size of fBuffer (at a given time)
   Int_t          fBufferLen;        ///< Current buffer length (<= fBufferSize)

   Long64_t       fBytesRead;        ///< Number of bytes read for this cache
   Long64_t       fBytesReadExtra;   ///< Number of extra bytes (overhead) read by the readahead buffer
   Int_t          fReadCalls;        ///< Number of read calls for this cache
   Long64_t       fNoCacheBytesRead; ///< Number of bytes read by basket to fill cached tree
   Int_t          fNoCacheReadCalls; ///< Number of read calls by basket to fill cached tree

   Bool_t         fAsyncReading;
   Bool_t         fEnablePrefetching;///< reading by prefetching asynchronously
   Bool_t         fEnableBuffering;  ///< buffer remote files on local disk

   Int_t          fNseek;            ///< Number of blocks to be prefetched
   Int_t          fNtot;             ///< Total size of prefetched blocks
   Int_t          fNb;               ///< Number of long buffers
   Int_t          fSeekSize;         ///< Allocated size of fSeek
   Long64_t      *fSeek;             ///<[fNseek] Position on file of buffers to be prefetched
   Long64_t      *fSeekSort;         ///<[fNseek] Position on file of buffers to be prefetched (sorted)
   Int_t         *fSeekIndex;        ///<[fNseek] sorted index table of fSeek
   Long64_t      *fPos;              ///<[fNb] start of long buffers
   Int_t         *fSeekLen;          ///<[fNseek] Length of buffers to be prefetched
   Int_t         *fSeekSortLen;      ///<[fNseek] Length of buffers to be prefetched (sorted)
   Int_t         *fSeekPos;          ///<[fNseek] Position of sorted blocks in fBuffer
   Int_t         *fLen;              ///<[fNb] Length of long buffers
   TFile         *fFile;             ///< Pointer to file
   char          *fBuffer;           ///<[fBufferSize] buffer of contiguous prefetched blocks
   Bool_t         fIsSorted;         ///< True if fSeek array is sorted
   Bool_t         fIsTransferred;    ///< True when fBuffer contains something valid
   Long64_t       fPrefetchedBlocks; ///< Number of blocks prefetched.

   //variables for the second block prefetched with the same semantics as for the first one
   Int_t          fBNseek;
   Int_t          fBNtot;
   Int_t          fBNb;
   Int_t          fBSeekSize;
   Long64_t      *fBSeek;        ///<[fBNseek]
   Long64_t      *fBSeekSort;    ///<[fBNseek]
   Int_t         *fBSeekIndex;   ///<[fBNseek]
   Long64_t      *fBPos;         ///<[fBNb]
   Int_t         *fBSeekLen;     ///<[fBNseek]
   Int_t         *fBSeekSortLen; ///<[fBNseek]
   Int_t         *fBSeekPos;     ///<[fBNseek]
   Int_t         *fBLen;         ///<[fBNb]
   Bool_t         fBIsSorted;
   Bool_t         fBIsTransferred;

   void SetEnablePrefetchingImpl(Bool_t setPrefetching = kFALSE); // Can not be virtual as it is called from the constructor.

private:
   TFileCacheRead(const TFileCacheRead &);            //cannot be copied
   TFileCacheRead& operator=(const TFileCacheRead &);

   std::unique_ptr<TFileBufferRead> fBufferFile;

public:
   TFileCacheRead();
   TFileCacheRead(TFile *file, Int_t buffersize, TObject *tree = 0);
   virtual ~TFileCacheRead();
   virtual Int_t       AddBranch(TBranch * /*b*/, Bool_t /*subbranches*/ = kFALSE) { return 0; }
   virtual Int_t       AddBranch(const char * /*branch*/, Bool_t /*subbranches*/ = kFALSE) { return 0; }
   virtual void        AddNoCacheBytesRead(Long64_t len) { fNoCacheBytesRead += len; }
   virtual void        AddNoCacheReadCalls(Int_t reads) { fNoCacheReadCalls += reads; }
   virtual void        Close(Option_t *option="");
   virtual Int_t       GetBufferSize() const { return fBufferSize; };
   virtual Long64_t    GetBytesRead() const { return fBytesRead; }
   virtual Long64_t    GetNoCacheBytesRead() const { return fNoCacheBytesRead; }
   virtual Long64_t    GetBytesReadExtra() const { return fBytesReadExtra; }
           TFile      *GetFile() const { return fFile; }   // Return the TFile being cached.
           Int_t       GetNseek() const { return fNseek; } // Return the number of blocks in the current cache.
           Int_t       GetNtot() const { return fNtot; }   // Return the total size of the prefetched blocks.
   virtual Int_t       GetReadCalls() const { return fReadCalls; }
   virtual Int_t       GetNoCacheReadCalls() const { return fNoCacheReadCalls; }
   virtual Int_t       GetUnzipBuffer(char ** /*buf*/, Long64_t /*pos*/, Int_t /*len*/, Bool_t * /*free*/) { return -1; }
           Long64_t    GetPrefetchedBlocks() const { return fPrefetchedBlocks; }
   virtual Bool_t      IsAsyncReading() const { return fAsyncReading; };
   virtual void        SetEnablePrefetching(Bool_t setPrefetching = kFALSE);
   virtual void        SetEnableBufferRead(bool setBuffer);
   virtual Bool_t      IsEnablePrefetching() const { return fEnablePrefetching; }
   virtual Bool_t      IsEnableBufferRead() const {return fBufferFile.get();}
   virtual Long_t      GetBufferedBlocks() const;
   virtual Bool_t      IsLearning() const {return kFALSE;}
   virtual Int_t       LearnBranch(TBranch * /*b*/, Bool_t /*subbranches*/ = kFALSE) { return 0; }
   virtual void        Prefetch(Long64_t pos, Int_t len);
   virtual void        Print(Option_t *option="") const;
   virtual Int_t       ReadBufferExt(char *buf, Long64_t pos, Int_t len, Int_t &loc);
   virtual Int_t       ReadBufferExtNormal(char *buf, Long64_t pos, Int_t len, Int_t &loc);
   virtual Int_t       ReadBufferExtPrefetch(char *buf, Long64_t pos, Int_t len, Int_t &loc);
   virtual Int_t       ReadBuffer(char *buf, Long64_t pos, Int_t len);
   virtual Int_t       SetBufferSize(Int_t buffersize);
   virtual void        SetFile(TFile *file, TFile::ECacheAction action = TFile::kDisconnect);
   virtual void        SetSkipZip(Bool_t /*skip*/ = kTRUE) {} // This function is only used by TTreeCacheUnzip (ignore it)
   virtual void        Sort();
   virtual void        SecondSort();                          //Method used to sort and merge the chunks in the second block
   virtual void        SecondPrefetch(Long64_t, Int_t);       //Used to add chunks to the second block
   virtual TFilePrefetch* GetPrefetchObj();
   virtual void        WaitFinishPrefetch();                  //Gracefully join the prefetching thread

   ClassDef(TFileCacheRead,2)  //TFile cache when reading
};

#endif
