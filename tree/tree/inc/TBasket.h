// @(#)root/tree:$Id$
// Author: Rene Brun   19/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBasket
#define ROOT_TBasket

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBasket                                                              //
//                                                                      //
// The TBasket objects are created at run time to collect TTree entries //
// in buffers. When a Basket is full, it is written to the file.        //
// The Basket is kept in memory if there is enough space.               //
//  (see the fMaxVirtualsize of TTree).                                 //
//                                                                      //
// The Basket class derives from TKey.                                  //
//////////////////////////////////////////////////////////////////////////


#include "TKey.h"

class TFile;
class TTree;
class TBranch;

class TBasket : public TKey {
friend class TBranch;

private:
   TBasket(const TBasket&);            ///< TBasket objects are not copiable.
   TBasket& operator=(const TBasket&); ///< TBasket objects are not copiable.

   // Internal corner cases for ReadBasketBuffers
   Int_t ReadBasketBuffersUnzip(char*, Int_t, Bool_t, TFile*);
   Int_t ReadBasketBuffersUncompressedCase();

   // Helper for managing the compressed buffer.
   void InitializeCompressedBuffer(Int_t len, TFile* file);

   // Handles special logic around deleting / reseting the entry offset pointer.
   void ResetEntryOffset();

   // Get entry offset as result of a calculation.
   Int_t *GetCalculatedEntryOffset();

   // Returns true if the underlying TLeaf can regenerate the entry offsets for us.
   Bool_t CanGenerateOffsetArray();

   // Manage buffer ownership.
   void   DisownBuffer();
   void   AdoptBuffer(TBuffer *user_buffer);

protected:
   Int_t       fBufferSize{0};                    ///< fBuffer length in bytes
   Int_t       fNevBufSize{0};                    ///< Length in Int_t of fEntryOffset OR fixed length of each entry if fEntryOffset is null!
   Int_t       fNevBuf{0};                        ///< Number of entries in basket
   Int_t       fLast{0};                          ///< Pointer to last used byte in basket
   Bool_t      fHeaderOnly{kFALSE};               ///< True when only the basket header must be read/written
   UChar_t     fIOBits{0};                        ///<!IO feature flags.  Serialized in custom portion of streamer to avoid forward compat issues unless needed.
   Bool_t      fOwnsCompressedBuffer{kFALSE};     ///<! Whether or not we own the compressed buffer.
   Bool_t      fReadEntryOffset{kFALSE};          ///<!Set to true if offset array was read from a file.
   Int_t      *fDisplacement{nullptr};            ///<![fNevBuf] Displacement of entries in fBuffer(TKey)
   Int_t      *fEntryOffset{nullptr};             ///<[fNevBuf] Offset of entries in fBuffer(TKey); generated at runtime.  Special value
                                                  /// of `-1` indicates that the offset generation MUST be performed on first read.
   TBranch    *fBranch{nullptr};                  ///<Pointer to the basket support branch
   TBuffer    *fCompressedBufferRef{nullptr};     ///<! Compressed buffer.
   Int_t       fLastWriteBufferSize[3] = {0,0,0}; ///<! Size of the buffer last three buffers we wrote it to disk
   Bool_t      fResetAllocation{false};           ///<! True if last reset re-allocated the memory
   UChar_t     fNextBufferSizeRecord{0};          ///<! Index into fLastWriteBufferSize of the last buffer written to disk
#ifdef R__TRACK_BASKET_ALLOC_TIME
   ULong64_t   fResetAllocationTime{0};           ///<! Time spent reallocating baskets in microseconds during last Reset operation.
#endif

   virtual void    ReadResetBuffer(Int_t basketnumber);

public:
   // The IO bits flag is to provide improved forward-compatibility detection.
   // Any new non-forward compatibility flags related serialization should be
   // added here.  When a new flag is added, set it in the kSupported field;
   //
   // The values and names of this (and EUnsupportedIOBits) enum need not be aligned
   // with the values of the various TIOFeatures enums, as there's a clean separation
   // between these two interfaces.  Practically, it is reasonable to keep them as aligned
   // as possible in order to avoid confusion.
   //
   // If (fIOBits & ~kSupported) is non-zero -- i.e., an unknown IO flag is set
   // in the fIOBits -- then the zombie flag will be set for this object.
   //
   enum class EIOBits : Char_t {
      // The following to bits are reserved for now; when supported, set
      // kSupported = kGenerateOffsetMap | kBasketClassMap
      kGenerateOffsetMap = BIT(0),
      // kBasketClassMap = BIT(1),
      kSupported = kGenerateOffsetMap
   };
   // This enum covers IOBits that are known to this ROOT release but
   // not supported; provides a mechanism for us to have experimental
   // changes that are not going go into a supported release.
   //
   // (kUnsupported | kSupported) should result in the '|' of all IOBits.
   enum class EUnsupportedIOBits : Char_t { kUnsupported = 0 };
   // The number of known, defined IOBits.
   static constexpr int kIOBitCount = 1;

   TBasket();
   TBasket(TDirectory *motherDir);
   TBasket(const char *name, const char *title, TBranch *branch);
   virtual ~TBasket();

   virtual void    AdjustSize(Int_t newsize);
   virtual void    DeleteEntryOffset();
   virtual Int_t   DropBuffers();
   TBranch        *GetBranch() const {return fBranch;}
           Int_t   GetBufferSize() const {return fBufferSize;}
           Int_t  *GetDisplacement() const {return fDisplacement;}
           Int_t *GetEntryOffset()
           {
              return R__likely(fEntryOffset != reinterpret_cast<Int_t *>(-1)) ? fEntryOffset : GetCalculatedEntryOffset();
           }
           Int_t   GetEntryPointer(Int_t Entry);
           Int_t   GetNevBuf() const {return fNevBuf;}
           Int_t   GetNevBufSize() const {return fNevBufSize;}
           Int_t   GetLast() const {return fLast;}
   virtual void    MoveEntries(Int_t dentries);
   virtual void    PrepareBasket(Long64_t /* entry */) {};
           Int_t   ReadBasketBuffers(Long64_t pos, Int_t len, TFile *file);
           Int_t   ReadBasketBytes(Long64_t pos, TFile *file);
   virtual void    WriteReset();

// Time spent reseting basket sizes (typically, at event cluster boundaries), in microseconds
#ifdef R__TRACK_BASKET_ALLOC_TIME
   ULong64_t       GetResetAllocationTime() const { return fResetAllocationTime; }
#endif
   // Count of resets performed of basket size.
   Bool_t          GetResetAllocationCount() const { return fResetAllocation; }

   Int_t           LoadBasketBuffers(Long64_t pos, Int_t len, TFile *file, TTree *tree = 0);
   Long64_t        CopyTo(TFile *to);

           void    SetBranch(TBranch *branch) { fBranch = branch; }
           void    SetNevBufSize(Int_t n) { fNevBufSize=n; }
   virtual void    SetReadMode();
   virtual void    SetWriteMode();
   inline  void    Update(Int_t newlast) { Update(newlast,newlast); };
   virtual void    Update(Int_t newlast, Int_t skipped);
   virtual Int_t   WriteBuffer();

   ClassDef(TBasket, 3); // the TBranch buffers
};

#endif
